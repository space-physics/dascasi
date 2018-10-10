#!/usr/bin/env python
"""
Reads DASC allsky cameras images in FITS formats into GeoData.
Run standalone from PlayDASC.py
"""
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
from typing import Union, Dict, Optional, Tuple, List
try:
    from skimage.transform import downscale_local_mean
except ImportError:
    downscale_local_mean = None

log = logging.getLogger('DASCutils-io')

# for NetCDF compression. too high slows down with little space savings.
ENC = {'zlib': True, 'complevel': 1, 'fletcher32': True}


def load(fin: Path,
         azelfn: Path = None,
         treq: np.ndarray = None,
         wavelenreq: List[str] = None,
         ofn: Path = None,
         verbose: bool = False) -> xarray.Dataset:
    """
    reads FITS images and spatial az/el calibration for allsky camera
    Bdecl is in degrees, from IGRF model

    odir: filename to write converted NetCDF4 (HDF5) video to.
    """
    flist = _slicereq(fin, treq, wavelenreq, verbose)

    time = []
    img: np.ndarray = []
    wavelen: np.ndarray = []

    warnings.filterwarnings('ignore', category=VerifyWarning)
    for fn in flist:
        try:
            im, t, w = _loadimg(fn)
        except (OSError, TypeError) as e:
            log.warning(f'{fn} has error {e}')
            continue

        img.append(im)
        time.append(t)
        wavelen.append(w)
    warnings.resetwarnings()

# %% collect output
    img = np.array(img)
    time = np.array(time).astype('datetime64[ns]')  # necessary for netcdf4 write
    wavelen = np.array(wavelen)

    data: xarray.Dataset = None
    for w in np.unique(wavelen):
        d = xarray.Dataset({w: (('time', 'y', 'x'), img[wavelen == w, ...])},
                           coords={'time': time[wavelen == w]})
        # 'y': range(img.shape[1]),
        # 'x': range(img.shape[2])})

        if data is None:
            data = d
        else:
            data = xarray.merge((data, d), join='outer')
# %% metadata
    data.attrs['filename'] = ' '.join((p.name for p in flist))
# %% camera location
    data = _camloc(flist[0], data)
# %% az / el
    data = _azel(azelfn, data)
# %% HDF5 write
    if isinstance(ofn, (str, Path)):
        ofn = Path(ofn).expanduser()
        if ofn.is_dir():
            ofn = ofn / f'{flist[0].name[:3]}_{data.time[0].values.astype(str)[:13]}.nc'
        else:
            ofn.parent.mkdir(exist_ok=True)

        print('writing', ofn)

        enc = {k: ENC for k in data.data_vars}
        data.to_netcdf(ofn, 'w', encoding=enc)

    return data


def _loadimg(fn: Path) -> Tuple[np.ndarray, datetime, str]:
    """
    DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
    just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
    since 14-bit {0,16384}.
    Further, there was a RAID failure that filled the data files with random values.
    Don Hampton says about 90% of data OK, but 10% NOK.
    """
    with fits.open(fn, mode='readonly', memmap=False) as h:
        assert h[0].header['BITPIX'] == 16, 'this function assumes unsigned 16-bit data'
        im = h[0].data.squeeze()  # Squeeze for old < 2011 files with 3-D, 1 image data.

    assert im.ndim == 2, 'one image at a time please'

    time = gettime(fn)
# %% rotation based on year
    k = -1 if time.year < 2012 else 1
    im = np.rot90(im, k=k)
# %% cleanup bad data
    im[im > 16384] = 0  # extreme, corrupted data
    im = im.clip(0, 16384).astype(np.uint16)  # discard bad values for 14-bit cameras.

    return im, time, getwavelength(fn)


def _slicereq(fin: Path, treq: np.ndarray,
              wavelenreq: List[str]=None, verbose: bool=False) -> List[Path]:

    forig = fin

    if isinstance(fin, (str, Path)):
        fin = Path(fin).expanduser()

    if isinstance(fin, Path):
        if fin.is_dir():
            flist = list(fin.glob('*.FITS'))
            if os.name != 'nt':
                flist += list(fin.glob('*.fits'))
            flist = sorted(flist)
        elif fin.is_file():
            flist = [fin]
        else:
            raise FileNotFoundError(forig)
# %% prefiltering files by user request for time or wavelength
    if treq is not None or wavelenreq is not None:
        time = []
        wavelen: np.ndarray = []
        flist1 = []
        for i, fn in enumerate(flist):
            try:
                time.append(gettime(fn))
                wavelen.append(getwavelength(fn))

                flist1.append(fn)
            except OSError:  # many corrupted files, accounted for by preallocated vectors
                pass

        time = np.atleast_1d(time)
        flist = np.atleast_1d(flist1)

        if len(wavelen) > 0:
            wavelen = np.atleast_1d(wavelen)
        else:  # old DASC files
            wavelen = None
# %% time request
    if treq is not None:
        if isinstance(treq, str):
            treq = parse(treq)
        elif isinstance(treq, (np.ndarray, tuple, list)):
            if isinstance(treq[0], str):
                treq = list(map(parse, treq))  # must have list()
            elif isinstance(treq[0], datetime):
                pass
            else:
                raise TypeError(f'not sure what time request you are making with {type(treq[0])}')

        treq = np.atleast_1d(treq)
# %% time slice
        if treq.size == 1:  # single frame
            i = abs(time-treq).argmin()  # nearest time, index number of file in flist
        elif treq.size == 2:  # frames within bounds
            i = (time >= treq[0]) & (time <= treq[1])  # boolean indexing

        if isinstance(i, np.ndarray) and not i.any():
            raise ValueError(f'no valid data found in {treq}')

        flist = np.atleast_1d(flist[i])
        wavelen = np.atleast_1d(wavelen[i])
        if len(flist) == 0:
            raise FileNotFoundError(f'no files found within time limits {treq}')
# %% wavelength slice
    if wavelenreq is not None and wavelen[0] != 'unknown':
        i = np.isin(wavelen, wavelenreq)
        flist = np.atleast_1d(flist[i])

        if len(flist) == 0:
            raise FileNotFoundError(f'no files found with wavelength(s) {wavelenreq}')
# %% iterate over image files
    if len(flist) == 0:
        raise FileNotFoundError(f'no DASC FITS files found in {forig}')

    if verbose:
        print('Number of files', len(flist), 'with wavelengths', np.unique(wavelen))

    return flist


def _camloc(fn: Path, data: xarray.Dataset) -> xarray.Dataset:
    """
    Camera altitude is not specified in the DASC files.
    This can be mitigated in the end user program with a WGS-84 height above
    ellipsoid lookup.
    """
    data.attrs['alt_m'] = np.nan

    lla = getcoords(fn)

    if lla is not None:
        data.attrs['lat'] = lla['lat']
        data.attrs['lon'] = lla['lon']
    else:
        data.attrs['lat'] = np.nan
        data.attrs['lon'] = np.nan

    return data


def _azel(azelfn: Optional[Path], data: xarray.Dataset) -> xarray.Dataset:

    if not azelfn:
        return data

    azel = loadcal(azelfn)

    wavelen = list(data.data_vars)

    imgshape = data[wavelen[0]].shape[1:]

    if azel['az'].shape != imgshape:
        downscale = (1, imgshape[0] // azel['az'].shape[0], imgshape[1] // azel['az'].shape[1])

        if downscale_local_mean is None:
            raise ImportError('pip install scikit-image')

        if downscale != 1:
            log.warning(f'downsizing images by factors of {downscale[1:]} to match calibration data')

        if len(wavelen) == 1 and wavelen[0] == 'unknown':
            if downscale != 1:
                data['unknown'] = (('time', 'y', 'x'), downscale_local_mean(data['unknown'], downscale))
            else:
                data['unknown'] = (('time', 'y', 'x'), data['unknown'])
        else:
            if downscale != 1:
                for w in np.unique(wavelen):
                    data[w] = downscale_local_mean(data[w], downscale)

    data.coords['az'] = (('y', 'x'), azel['az'])
    data.coords['el'] = (('y', 'x'), azel['el'])

    return data


def loadcal(azelfn: Path) -> xarray.Dataset:
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
        raise ValueError('Az and El are the same file!')

    with fits.open(azfn, mode='readonly') as h:
        az = h[0].data
    bad = az == 0

    with fits.open(elfn, mode='readonly') as h:
        el = h[0].data
    bad &= el == 0.

    el[bad] = np.nan
    az[bad] = np.nan

    assert np.nanmax(el) <= 90 and np.nanmin(el) >= 0, '0 < elevation < 90 degrees.'
    assert np.nanmax(az) <= 360 and np.nanmin(az) >= 0, '0 < azimuth < 360 degrees.'

    azel = xarray.Dataset({'el': (('y', 'x'), el),
                           'az': (('y', 'x'), az)})

    return azel


def gettime(fn: Path) -> datetime:
    """ returns time of DASC frame in file (assumes one frame per file)"""
    with fits.open(fn, mode='readonly') as h:
        try:
            t = parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART'])
            #   expstart = parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART'])
#               time.append((expstart, expstart + timedelta(seconds=h[0].header['EXPTIME']))) #EXPTIME is in seconds
        except KeyError:
            t = parse(h[0].header['FRAME'])

    return t


def getwavelength(fn: Path) -> str:
    """ returns optical wavelength [nm] of DASC frame in file (assumes one frame per file)"""

    with fits.open(fn) as h:
        try:
            w = h[0].header['FILTWAV']
        except KeyError:
            w = 'unknown'

    return w


def getcoords(fn: Path) -> Optional[Dict[str, float]]:
    """ get lat, lon from DASC header"""

    latlon: Union[None, Dict[str, float]]

    with fits.open(fn) as h:
        try:
            latlon = {'lat': h[0].header['GLAT'],
                      'lon': h[0].header['GLON']}
        except KeyError:
            if fn.name[:3] == 'PKR':
                latlon = {'lat': 65.126,
                          'lon': -147.479}
            else:
                latlon = None

    return latlon


def stem2fn(stem: Path) -> Tuple[Path, Path]:
    """if user specifies the stem to Az,El, generate the az, el filenames"""
    assert isinstance(stem, (str, Path))

    stem = Path(stem).expanduser()

    if not stem.parent.is_dir():
        raise FileNotFoundError(f'{stem.parent} is not a directory')
    elif stem.is_file():
        raise OSError(f'need to specify stem, not whole filename {stem}')

    azfn = stem.parent / (stem.name + '_AZ_10deg.fits')
    elfn = stem.parent / (stem.name + '_EL_10deg.fits')

    if not azfn.is_file() or not elfn.is_file():
        raise FileNotFoundError(f'did not find {azfn} \n {elfn}')

    return azfn, elfn
