#!/usr/bin/env python
"""
Reads DASC allsky cameras images in FITS formats into GeoData.
Run standalone from PlayDASC.py
"""
import sys
from pathlib import Path
import warnings  # corrupt FITS files let off a flood of AstroPy warnings
from astropy.io.fits.verify import VerifyWarning
import logging
from astropy.io import fits
import numpy as np
from datetime import datetime
from dateutil.parser import parse
import xarray
from typing import Union, Sequence, Dict, Optional, Tuple
try:
    from skimage.transform import downscale_local_mean
except ImportError:
    downscale_local_mean = None
    
from .processing import interpolateCoordinate, circular2lla, datetime2posix

log = logging.getLogger('DASCutils-io')

def load(fin: Union[Path, Sequence[Path]],
         azelfn: Union[Path, Sequence[Path]]=None,
         treq: np.ndarray=None,
         wavelenreq: list=None, verbose: bool=False,
         coordinate: str='polar',
         mapping_altitude: Union[int,float]=None,
         ofn: Union[Path,str]=None) -> xarray.Dataset:
    """
    reads FITS images and spatial az/el calibration for allsky camera
    Bdecl is in degrees, from IGRF model
    """
    forig = fin
    warnings.filterwarnings('ignore', category=VerifyWarning)

    if isinstance(fin, (str, Path)):
        fin = Path(fin).expanduser()

    if isinstance(fin, Path):
        if fin.is_dir():
            flist = list(fin.glob('*.FITS'))
            if sys.platform != 'win32':
                flist += list(fin.glob('*.fits'))
            flist = sorted(flist)
        elif fin.is_file():
            flist = [fin]
        else:
            raise FileNotFoundError(f'not sure what {flist} is')
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
    if wavelenreq is not None and wavelen is not None:
        i = np.isin(wavelen, wavelenreq)
        flist = flist[i]

        if len(flist) == 0:
            raise FileNotFoundError(f'no files found with wavelength(s) {wavelenreq}')
# %% iterate over image files
    if len(flist) == 0:
        raise FileNotFoundError(f'no DASC FITS files found in {forig}')

    if verbose:
        print('Number of files', len(flist), 'with wavelengths', np.unique(wavelen))

    time = []
    img: np.ndarray = []
    wavelen = []
    lla = None
    # Get a rotation angle, different for different cameras
    k = -1 if (int(fn.name[14:18]) < 2012) else 1
    for i, fn in enumerate(flist):
        try:
            if i == 0:
                lla = getcoords(fn)

            with fits.open(fn, mode='readonly', memmap=False) as h:
                if verbose and i == 0:
                    print(h[0].header)

                assert h[0].header['BITPIX'] == 16, 'this function assumes unsigned 16-bit data'

                """
                DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
                just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
                since 14-bit \in {0,16384}.
                Further, there was a RAID failure that filled the data files with random values.
                Don Hampton says about 90% of data OK, but 10% NOK.
                """

                im = h[0].data.squeeze()  # Squeeze for old < 2011 files with 3-D, 1 image data.
                assert im.ndim == 2, 'one image at a time please'
                
                im = np.rot90(im, k=k)  # NOTE: rotate by -1 to match online AVIs from UAF website.
                im[im > 16384] = 0  # extreme, corrupted data
                im = im.clip(0, 16384).astype(np.uint16)  # discard bad values for 14-bit cameras.

            img.append(im)

            # read headers at the end, as failures virtually always happens when reading the image (truncated files?) vs. header
            time.append(gettime(fn))
            wavelen.append(getwavelength(fn))

        except (OSError, TypeError) as e:
            log.warning(f'{fn} has error {e}')

# %% collect output
    img = np.array(img)

    time = np.array(time)  # this prevents xarray from using nanaseconds M8 datetime64 that is annoying.
    if len(wavelen) == 0 or wavelen[0] is None:
        wavelen = None

    if wavelen is None:
        data = xarray.Dataset({'000': (('time', 'y', 'x'), img)},
                              coords={'time': time})
        if data.time.size > 1:  # more than 1 image
            data.attrs['cadence'] = round( ((time[1]-time[0]).total_seconds()),2)  # NOTE: assumes uniform kinetic rate
    else:
        data = None
        cadence = {}
        for w in np.unique(wavelen):
            d = xarray.Dataset({str(w): (('time', 'y', 'x'), img[wavelen == w, ...])},
                               coords={'time': time[wavelen == w]})
            if d.time.size > 1:  # more than 1 image
                cadence[w] = round(((d.time[1] - d.time[0]).item()*1e-9),2)  # NOTE: assumes uniform kinetic rate
            else:
                cadence[w] = 0 #

            if data is None:
                data = d
            else:
                data = xarray.merge((data, d), join='outer')

        data.attrs['cadence'] = cadence
# %% coordinates of camera
    """
    Camera altitude is not specified in the DASC files.
    This can be mitigated in the end user program with a WGS-84 height above
    ellipsoid lookup.
    """
    data.attrs['alt0'] = 0
    
    if lla is not None:
        data.attrs['lat0'] = lla['lat']
        data.attrs['lon0'] = lla['lon']
    else:
        data.attrs['lat0'] = None
        data.attrs['lon0'] = None
# %% az / el
    if azelfn is not None:
        azel = loadcal(azelfn)
        if azel['az'].shape != im.shape:
            downscale = (1, im.shape[0] // azel['az'].shape[0], im.shape[1] // azel['az'].shape[1])

            if downscale_local_mean is None:
                raise ImportError('pip install scikit-image')

            if downscale != 1:
                log.warning(f'downsizing images by factors of {downscale[1:]} to match calibration data')

            if wavelen is None:
                if downscale != 1:
                    data['000'] = (('time', 'y', 'x'), downscale_local_mean(data['000'], downscale))
                else:
                    data['000'] = (('time', 'y', 'x'), data['000'])
            else:
                if downscale != 1:
                    for w in np.unique(wavelen):
                        data[w] = downscale_local_mean(data[w], downscale)

        data['az'] = azel['az']
        data['el'] = azel['el']
        
# %% az/el -> lat/lon
    if coordinate == 'wsg':
        if mapping_altitude is None:
            mapping_altitude = 100
            print ('Attribute mapping altitude was not set! Deafult is set to 100 km.')

        # Get rid of NaNs in the coordinates' arrays
        eli = interpolateCoordinate(azel['el'].values, N = azel['el'].values.shape[0], method = 'nearest')
        azi = interpolateCoordinate(azel['az'].values, N = azel['el'].values.shape[0], method = 'nearest')
        # Convert Coordinates to WSG84
        lat, lon, alt = circular2lla(az = azi, el = eli, lat0 = lla['lat'], 
                                     lon0 = lla['lon'], alt0 = 0,
                                     mapping_altitude = mapping_altitude)
        data.coords['lat'] = (('x','y'), lat)
        data.coords['lon'] = (('y','x'), lon)
        data.attrs['alt_m'] = mapping_altitude
    data.attrs['filename'] = ' '.join((p.name for p in flist))
    data.attrs['wavelength'] = wavelen if isinstance(wavelen, (int,float)) else '000'

# %% Save to netCDF?
    if ofn is None:
        return data
    else:
        # Convert datetime to posix. netCDf cannot save datetime
        posix_time = datetime2posix(data.time.values)
        data['time'] = posix_time
        if wavelen is not None:
            data.attrs['cadence'] = cadence[wavelen[0]]
        # To HDF
        try: 
            save2HDF(data=data, fn_out=ofn)
            return data
        except Exception as e:
            return data
            raise (e)


def loadcal(azelfn: Union[Path, Sequence[Path]]) -> xarray.Dataset:
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


def getwavelength(fn: Path) -> Optional[int]:
    """ returns optical wavelength [nm] of DASC frame in file (assumes one frame per file)"""

    w: Union[None, int]

    with fits.open(fn) as h:
        try:
            w = int(h[0].header['FILTWAV'])
        except KeyError:
            w = None

    return w


def getcoords(fn: Path) -> Optional[Dict[str, float]]:
    """ get lat, lon from DASC header"""

    latlon: Union[None, Dict[str, float]]

    with fits.open(fn) as h:
        try:
            latlon = {'lat': h[0].header['GLAT'],
                      'lon': h[0].header['GLON']}
        except KeyError:
            if 'PKR' in fn.name:
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

def save2HDF(data: xarray.Dataset=None,
            fn_out: Union[str, Path]=None):
    ENC = {'zlib': True, 'complevel': 1, 'fletcher32': True}
    if data is None:
        raise ('Input data invalid')
    if (fn_out is '') or (fn_out is None):
        print ('Path/output file not given. Default is desktop!')
        fn_out = Path.home() / 'Desktop' / 'dasc_out.nc'
    
    if isinstance(fn_out, str): fn_out = Path(fn_out).expanduser()
    
    if fn_out.is_dir(): fn_out = fn_out / 'dasc_out.nc'
    
    if (fn_out.suffix != '.nc') or (fn_out.suffix != '.h5') or (fn_out.suffix != '.hdf5'):
        fn_out = fn_out.with_suffix('.nc')
        
    # Write a timestamp of conversion
    data.attrs['converted'] = str(datetime.now())
    
    enc = {k:ENC for k in data.data_vars}
    print ('Saving data into: ', fn_out)
    data.to_netcdf(fn_out, mode='w', encoding=enc, group='DASC')