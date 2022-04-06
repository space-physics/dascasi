"""
functions using HDF5
"""

from __future__ import annotations
import h5py
import typing as T
from pathlib import Path
import numpy as np
from datetime import datetime
import xarray

from .utils import get_time_slice


def save_hdf5(imgs: dict[str, T.Any], outfile: Path):
    print("writing image stack to", outfile)

    with h5py.File(outfile, "w") as f:
        f["wavelengths"] = imgs["wavelengths"].astype(np.string_)
        for p in ["alt0", "lat0", "lon0", "az", "el"]:
            if p in imgs:
                f[f"/camera/{p}"] = imgs[p]
        for wl in imgs["wavelengths"]:
            h = f.create_dataset(
                f"/{wl}/imgs",
                data=imgs[wl],
                compression="gzip",
                compression_opts=1,
                # because of high entropy data,
                # higher compression didn't significantly shrink file size
                chunks=(1, *imgs[wl].shape[1:]),
                shuffle=True,
                fletcher32=True,
            )
            # metadata to show this is an image stack
            h.attrs["CLASS"] = np.string_("IMAGE")
            h.attrs["IMAGE_VERSION"] = np.string_("1.2")
            h.attrs["IMAGE_SUBCLASS"] = np.string_("IMAGE_GRAYSCALE")
            h.attrs["DISPLAY_ORIGIN"] = np.string_("LL")
            h.attrs["IMAGE_WHITE_IS_ZERO"] = np.uint8(0)
            # time vector
            f[f"/{wl}/time"] = imgs[wl].time.values.astype(np.string_)
            # camera location
            if "lat" in imgs[wl].coords:
                f[f"/{wl}/lat"] = imgs[wl].lat
                f[f"/{wl}/lon"] = imgs[wl].lon


def load_hdf5(
    filename: Path, treq: list[datetime] = None, wavelenreq: list[str] = None
) -> dict[str, T.Any]:

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
                coords={
                    "time": time[i],
                    "y": range(h[f"/{wl}/imgs"].shape[1]),
                    "x": range(h[f"/{wl}/imgs"].shape[2]),
                },
                dims=["time", "y", "x"],
            )

    return imgs
