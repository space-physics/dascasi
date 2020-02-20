"""
functions using HDF5
"""
import h5py
import typing as T
from pathlib import Path
import numpy as np
from datetime import datetime
import xarray

from .utils import get_time_slice


def save_hdf5(imgs: T.Dict[str, T.Any], outfile: Path):
    print("writing image stack to", outfile)

    with h5py.File(outfile, "w") as h:
        h["wavelengths"] = imgs["wavelengths"].astype(np.string_)
        for p in ["alt0", "lat0", "lon0", "az", "el"]:
            if p in imgs:
                h[f"/camera/{p}"] = imgs[p]
        for wl in imgs["wavelengths"]:
            h.create_dataset(
                f"/{wl}/imgs",
                data=imgs[wl],
                compression="gzip",
                compression_opts=1,
                chunks=(1, *imgs[wl].shape[1:]),
                shuffle=True,
                fletcher32=True,
            )
            h[f"/{wl}/time"] = imgs[wl].time.values.astype(np.string_)
            if "lat" in imgs[wl].coords:
                h[f"/{wl}/lat"] = imgs[wl].lat
                h[f"/{wl}/lon"] = imgs[wl].lon


def load_hdf5(filename: Path, treq: T.Sequence[datetime] = None, wavelenreq: T.Sequence[str] = None) -> T.Dict[str, T.Any]:

    imgs = {}

    with h5py.File(filename, "r") as h:
        for p in ["alt0", "lat0", "lon0"]:
            imgs[p] = h[f"camera/{p}"][()]
        if "az" in h:
            imgs["az"] = h["camera/az"][:]
            imgs["el"] = h["camera/el"][:]

        wavelen = h["wavelengths"][:].astype(str) if wavelenreq is None else wavelenreq
        imgs["wavelengths"] = wavelen

        for wl in wavelen:
            time = np.asarray(h[f"/{wl}/time"]).astype("datetime64[us]")
            i = get_time_slice(time, treq)
            imgs[wl] = xarray.DataArray(
                data=h[f"/{wl}/imgs"][i, ...],
                name=wl,
                coords={"time": time[i], "y": range(h[f"/{wl}/imgs"].shape[1]), "x": range(h[f"/{wl}/imgs"].shape[2])},
                dims=["time", "y", "x"],
            )

    return imgs
