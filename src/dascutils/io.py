#!/usr/bin/env python
"""
Reads DASC allsky cameras images in FITS formats into GeoData.
Run standalone from PlayDASC.py
"""

from __future__ import annotations
import os
from pathlib import Path
import warnings  # corrupt FITS files let off a flood of AstroPy warnings
from astropy.io.fits.verify import VerifyWarning
import logging
from astropy.io import fits
import numpy as np
from datetime import datetime
from dateutil.parser import parse
import xarray
import typing as T
import pymap3d as pm

from .projection import interpolateCoordinate, interpSpeedUp
from .hdf5 import load_hdf5
from .utils import get_time_slice

try:
    from skimage.transform import downscale_local_mean
except ImportError:
    downscale_local_mean = None


log = logging.getLogger("DASCutils-io")


def load(
    fin: Path,
    azelfn: Path = None,
    treq: list[datetime] = None,
    wavelenreq: list[str] = None,
    wavelength_altitude_km: dict[str, float] = None,
) -> dict[str, T.Any]:
    """
    reads FITS images and spatial az/el calibration for allsky camera
    Bdecl is in degrees, from IGRF model
    """
    fin = Path(fin).expanduser()
    if fin.is_file() and fin.suffix in (".h5", ".hdf5"):
        return load_hdf5(fin, treq, wavelenreq)

    flist = _slicereq(fin, treq, wavelenreq)
    if not flist:
        raise FileNotFoundError(f"No files found in {fin}")

    # %% load data from good files, discarding bad
    imgs = _sift(flist)
    # %% camera location
    imgs = _camloc(imgs, flist[0].parent)
    # %% az / el
    imgs = _azel(azelfn, imgs)
    # %% projections
    imgs = _project(imgs, wavelength_altitude_km)

    return imgs


def _sift(flist: list[Path]) -> dict[str, xarray.DataArray]:
    """find good files from bad files and preload the data"""

    time = []
    files = []  # to only keep filesnames for "good" files
    img = []
    wavelen: list[str] = []

    warnings.filterwarnings("ignore", category=VerifyWarning)
    for fn in flist:
        try:
            im, t, w = _loadimg(fn)
        except (OSError, TypeError, ValueError) as e:
            log.warning(f"{fn} has error {e}")
            continue

        files.append(fn)
        img.append(im)
        time.append(t)
        wavelen.append(w)
    warnings.resetwarnings()

    return _collect(files, img, time, wavelen)


def _project(
    imgs: dict[str, xarray.DataArray], wavelength_altitude_km: dict[str, float]
) -> dict[str, xarray.DataArray]:
    """project image stack to specified per-wavelength altitudes"""

    if wavelength_altitude_km is None:
        return imgs

    eli = interpolateCoordinate(imgs["el"], method="nearest")
    azi = interpolateCoordinate(imgs["az"], method="nearest")

    for wl, mapalt_km in wavelength_altitude_km.items():
        if wl not in imgs:
            continue
        lat, lon, _ = pm.aer2geodetic(
            azi,
            eli,
            mapalt_km * 1e3 / np.sin(np.radians(eli)),
            imgs["lat0"],
            imgs["lon0"],
            imgs["alt0"],
        )

        mapped_lon, mapped_lat, mapped_img = interpSpeedUp(
            x_in=lon, y_in=lat, image=imgs[wl].values
        )
        imgs[wl].data = mapped_img

        # lat, lon cannot be dimensions here because they're each dynamic in 2-D
        imgs[wl].coords["lat"] = (("y", "x"), mapped_lat)
        imgs[wl].coords["lon"] = (("y", "x"), mapped_lon)
        imgs[wl].attrs["mapping_alt_km"] = mapalt_km

    return imgs


def _collect(
    files: list[Path], img: list, time: list[datetime], wavelen: list[str]
) -> dict[str, xarray.DataArray]:
    """assemble image stack into dict of xarray.DataArray"""
    img = np.array(img)
    time = np.array(time)
    wavelen = np.array(wavelen)
    files = np.asarray(files)  # for boolean indexing

    imgs = {"wavelengths": np.unique(wavelen)}
    for w in imgs["wavelengths"]:
        i = wavelen == w
        imgs[w] = xarray.DataArray(
            data=img[i, ...],
            name=w,
            coords={"time": time[i], "y": range(img.shape[1]), "x": range(img.shape[2])},
            dims=["time", "y", "x"],
        )
        imgs[w].attrs["filename"] = [p.name for p in files[i]]

    return imgs


def _loadimg(fn: Path) -> tuple[T.Any, datetime, str]:
    """
    DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
    just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
    since 14-bit {0,16384}.
    Further, there was a RAID failure that filled the data files with random values.
    Don Hampton says about 90% of data OK, but 10% NOK.
    """
    with fits.open(fn, mode="readonly", memmap=False) as h:
        assert h[0].header["BITPIX"] == 16, "this function assumes unsigned 16-bit data"
        im = h[0].data.squeeze()  # Squeeze for old < 2011 files with 3-D, 1 image data.

    assert im.ndim == 2, "one image at a time please"

    time = gettime(fn)
    # %% rotation based on year
    k = -1 if time.year < 2012 else 1
    im = np.rot90(im, k=k)
    # %% cleanup bad data
    im[im > 16384] = 0  # extreme, corrupted data
    im = im.clip(0, 16384).astype(np.uint16)  # discard bad values for 14-bit cameras.

    return im, time, getwavelength(fn)


def _slicereq(fin: Path, treq: list[datetime], wavelenreq: list[str] = None) -> list[Path]:
    """given user parameters, determine slice for image stack vs. wavelength and time"""

    if fin.is_dir():
        flist = list(fin.glob("*.FITS"))
        if os.name != "nt":
            flist += list(fin.glob("*.fits"))
    elif fin.is_file():
        flist = [fin]
    else:
        raise FileNotFoundError(f"{fin} is not a file or directory")

    # we may filter by time, so put in time order
    flist = sorted(flist)

    if treq is None and wavelenreq is None:
        return flist
    # %% prefiltering files by user request for time or wavelength
    time = []
    wavelen = []
    flist1 = []
    for fn in flist:
        try:
            time.append(gettime(fn))
            wavelen.append(getwavelength(fn))
            flist1.append(fn)
        except OSError:
            # there are many corrupted files
            pass

    flist = flist1
    # %% time request
    i = get_time_slice(time, treq)
    flist = flist[i]
    wavelen = wavelen[i]
    # %% wavelength slice
    if wavelenreq is not None:
        j = np.isin(wavelen, wavelenreq)
        flist = np.asarray(flist)[j].tolist()

    return flist


def _camloc(imgs: dict[str, T.Any], path: Path) -> dict[str, T.Any]:
    """
    Camera altitude is not specified in the DASC files.
    This can be mitigated in the end user program with a WGS-84 height above
    ellipsoid lookup.
    """

    imgs["alt0"] = 0.0
    # arbitrarily pick a filename as camera is stationary
    filename = path / imgs[imgs["wavelengths"][0]].filename[0]
    imgs.update(getcoords(filename))

    return imgs


def _azel(azelfn: Path, data: dict[str, T.Any]) -> dict[str, T.Any]:

    if not azelfn:
        return data

    azel = loadcal(azelfn)

    wavelen = data["wavelengths"]

    imgshape = data[wavelen[0]].shape[1:]

    if azel["az"].shape != imgshape:
        downscale = (1, imgshape[0] // azel["az"].shape[0], imgshape[1] // azel["az"].shape[1])

        if downscale_local_mean is None:
            raise ImportError("pip install scikit-image")

        if downscale != 1:
            log.warning(
                f"downsizing images by factors of {downscale[1:]} to match calibration data"
            )

        if len(wavelen) == 1 and wavelen[0] == "0000":
            if downscale != 1:
                data["0000"] = downscale_local_mean(data["0000"], downscale)
        else:
            if downscale != 1:
                for w in np.unique(wavelen):
                    data[w] = downscale_local_mean(data[w], downscale)

    data.update(azel)

    return data


def loadcal(azelfn: Path) -> dict[str, T.Any]:
    """Load DASC plate scale (degrees/pixel)"""
    if isinstance(azelfn, (str, Path)):
        azfn, elfn = stem2fn(azelfn)
    elif len(azelfn) == 1:
        azfn, elfn = stem2fn(azelfn[0])
    else:
        azfn, elfn = azelfn

    azfn = Path(azfn).expanduser()
    elfn = Path(elfn).expanduser()

    if azfn.samefile(elfn):
        raise ValueError("Az and El are the same file!")

    with fits.open(azfn, mode="readonly") as h:
        az = h[0].data
    bad = az == 0

    with fits.open(elfn, mode="readonly") as h:
        el = h[0].data
    bad &= el == 0.0

    el[bad] = np.nan
    az[bad] = np.nan

    assert np.nanmax(el) <= 90 and np.nanmin(el) >= 0, "0 < elevation < 90 degrees."
    assert np.nanmax(az) <= 360 and np.nanmin(az) >= 0, "0 < azimuth < 360 degrees."

    azel = {"el": el, "az": az}

    return azel


def gettime(fn: Path) -> datetime:
    """returns time of DASC frame in file (assumes one frame per file)"""
    with fits.open(fn, mode="readonly") as h:
        try:
            t = parse(h[0].header["OBSDATE"] + "T" + h[0].header["OBSSTART"])
            #   expstart = parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART'])
        #               time.append((expstart, expstart + timedelta(seconds=h[0].header['EXPTIME']))) #EXPTIME is in seconds
        except KeyError:
            t = parse(h[0].header["FRAME"])

    return t


def getwavelength(fn: Path) -> str:
    """returns optical wavelength [nm] of DASC frame in file.
    Assumes one frame per file."""

    with fits.open(fn) as h:
        try:
            w = h[0].header["FILTWAV"]
        except KeyError:
            # extract from filename (older DASC FITS files have this convention too)
            w = fn.stem[9:13]
            try:
                int(w)
            except ValueError:
                w = "0000"

    return w


def getcoords(fn: Path) -> dict[str, float]:
    """get lat, lon from DASC header"""

    with fits.open(fn) as h:
        try:
            latlon = {"lat0": h[0].header["GLAT"], "lon0": h[0].header["GLON"]}
        except KeyError:
            if fn.name[:3] == "PKR":
                latlon = {"lat0": 65.126, "lon0": -147.479}
            else:
                latlon = {"lat0": np.nan, "lon0": np.nan}

    return latlon


def stem2fn(stem: Path) -> tuple[Path, Path]:
    """if user specifies the stem to Az,El, generate the az, el filenames"""
    assert isinstance(stem, (str, Path))

    stem = Path(stem).expanduser()

    if not stem.parent.is_dir():
        raise FileNotFoundError(f"{stem.parent} is not a directory")
    elif stem.is_file():
        raise OSError(f"need to specify stem, not whole filename {stem}")

    stemname = stem.name.rstrip("_")

    azfn = stem.parent / (stemname + "_Az.FIT")
    elfn = stem.parent / (stemname + "_El.FIT")

    if not azfn.is_file() or not elfn.is_file():
        raise FileNotFoundError(f"did not find {azfn} \n {elfn}")

    return azfn, elfn
