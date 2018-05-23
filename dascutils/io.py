#!/usr/bin/env python
"""
Reads DASC allsky cameras images in FITS formats into GeoData.
Run standalone from PlayDASC.py
"""
from pathlib import Path
import warnings # corrupt FITS files let off a flood of AstroPy warnings
from astropy.io.fits.verify import VerifyWarning
import logging
from astropy.io import fits
import numpy as np
from datetime import datetime
from dateutil.parser import parse
import xarray
from typing import Tuple


def load(flist:list, azfn:Path=None, elfn:Path=None, treq:list=None,
         wavelenreq:list=None, verbose:bool=False) -> xarray.Dataset:
    """
    reads FITS images and spatial az/el calibration for allsky camera
    Bdecl is in degrees, from IGRF model
    """
    warnings.filterwarnings('ignore', category=VerifyWarning)

    if isinstance(flist,(str,Path)):
        flist = Path(flist).expanduser()

    if isinstance(flist,Path):
        if flist.is_dir():
            flist = list(flist.glob('*.FITS')) + list(flist.glob('*.fits'))
            flist = sorted(flist)
        elif flist.is_file():
            flist = [flist]
        else:
            raise FileNotFoundError(f'not sure what {flist} is')
# %% prefiltering files by user request for time or wavelength
    if treq is not None or wavelenreq is not None:
        time = []; wavelen=[]; flist1 = []
        for i,fn in enumerate(flist):
            try:
                with fits.open(fn, mode='readonly') as h:
                    try:
                        time.append(parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART']))
                    except KeyError:
                        time.append(parse(h[0].header['FRAME']))

                    try:
                        wavelen.append(int(h[0].header['FILTWAV']))
                    except KeyError:
                        pass

                flist1.append(fn)
            except OSError: #many corrupted files, accounted for by preallocated vectors
                pass

        time = np.atleast_1d(time)
        flist = np.atleast_1d(flist1)

        if len(wavelen) > 0:
            wavelen = np.atleast_1d(wavelen)
        else: # old DASC files
            wavelen = None
#%% time request
    if treq is not None:
        if isinstance(treq,str):
            treq = parse(treq)
        elif isinstance(treq,(np.ndarray,tuple,list)):
            if isinstance(treq[0],str):
                treq = list(map(parse,treq)) # must have list()
            elif isinstance(treq[0],datetime):
                pass
            else:
                raise TypeError(f'not sure what time request you are making with {type(treq[0])}')

        treq = np.atleast_1d(treq)
# %% time slice
        if treq.size == 1: # single frame
            i = abs(time-treq).argmin() # nearest time, index number of file in flist
        elif treq.size == 2: #frames within bounds
            i = (time>=treq[0]) & (time<=treq[1]) # boolean indexing

        if isinstance(i,np.ndarray) and not i.any():
            raise ValueError(f'no valid data found in {treq}')

        flist   = np.atleast_1d(flist[i])
        wavelen = np.atleast_1d(wavelen[i])
        if len(flist)==0:
            raise FileNotFoundError(f'no files found within time limits {treq}')
# %% wavelength slice
    if wavelenreq is not None and wavelen is not None:
        i = np.isin(wavelen, wavelenreq)
        flist = np.atleast_1d(flist[i])

        if len(flist)==0:
            raise FileNotFoundError(f'no files found with wavelength(s) {wavelenreq}')
#%% iterate over image files
    if verbose:
        print('Number of files',len(flist),'with wavelengths',np.unique(wavelen))

    time = []; img= [];  wavelen = []; lla=None
    for i,fn in enumerate(flist):
        try:
            with fits.open(fn, mode='readonly', memmap=False) as h:
                if verbose and i==0:
                    print(h[0].header)

                assert h[0].header['BITPIX']==16,'this function assumes unsigned 16-bit data'
                try:
                     time.append(parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART']))
#                    expstart = parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART'])
#                    time.append((expstart, expstart + timedelta(seconds=h[0].header['EXPTIME']))) #EXPTIME is in seconds
                except KeyError: #old DASC files < 2011
                    time.append(parse(h[0].header['FRAME']))

                try:
                    wavelen.append(int(h[0].header['FILTWAV']))
                except KeyError:
                    pass

                try:
                    lla={'lat':h[0].header['GLAT'],
                         'lon':h[0].header['GLON']}
                except KeyError:
                    pass
                """
                DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
                just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
                since 14-bit \in {0,16384}.
                These files do not have a BZERO value. Someday when they're written correctly this
                code may need updating.
                Further, there was a RAID failure that filled the data files with random values.
                Don Hampton says about 90% of data OK, but 10% NOK.
                """
                I = h[0].data.squeeze()  # Squeeze for old < 2011 files with 3-D, 1 image data.
                assert I.ndim == 2,'one image at a time please'

                I = np.rot90(I, k=1) # NOTE: rotate by -1 to match online AVIs from UAF website.
         #       if not 'BZERO' in h[0].header:
         #           I[I>16384] = 0 #extreme, corrupted data
         #           I = I.clip(0,16384).astype(np.uint16) #discard bad values for 14-bit cameras.

            img.append(I)

        except (OSError,TypeError) as e:
            logging.warning(f'{fn} has error {e}')

# %% collect output
    img = np.array(img)
    time = np.array(time) # this prevents xarray from using nanaseconds M8 datetime64 that is annoying.
    if len(wavelen) == 0:
        wavelen = None


    if wavelen is None:
        data = xarray.Dataset({'white': (('time','y','x'), img)},
                               coords={'time':time})
    else:
        data = None
        for w in np.unique(wavelen):
            d = xarray.Dataset({w:(('time','y','x'),img[wavelen==w,...])},
                               coords={'time':time[wavelen==w]})
            if data is None:
                data = d
            else:
                data = xarray.merge((data,d))


 #                         attrs={#'timeend':time[:,1],
    if lla is not None:
        data.attrs['lat']=lla['lat']
        data.attrs['lon']=lla['lon']

    if azfn is not None and elfn is not None:
        az,el = loadcal(azfn, elfn)
        data['az'] = az
        data['el'] = el

    data.attrs['filename'] = ' '.join((p.name for p in flist))
    data.attrs['wavelength'] = wavelen

    return data


def loadcal(azfn:Path, elfn:Path) -> Tuple[np.ndarray,np.ndarray]:
    """Load DASC plate scale (degrees/pixel)"""

    azfn = Path(azfn).expanduser()
    elfn = Path(elfn).expanduser()

    if azfn.samefile(elfn):
        raise ValueError('Az and El are the same file!')

    with fits.open(azfn, mode='readonly') as h:
        az = h[0].data
    bad = az==0

    with fits.open(elfn, mode='readonly') as h:
        el = h[0].data
    bad &= el==0.

    el[bad] = np.nan
    az[bad] = np.nan

    assert np.nanmax(el) <= 90  and np.nanmin(el) >= 0, '0 < elevation < 90 degrees.'
    assert np.nanmax(az) <= 360 and np.nanmin(az) >= 0, '0 < azimuth < 360 degrees.'

    el = (('y','x'), el)
    az = (('y','x'), az)

    return az,el
