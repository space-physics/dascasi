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
from datetime import timedelta
from dateutil.parser import parse
import xarray


def load(flist:list, azfn:Path=None, elfn:Path=None, treq:list=None, wavelenreq:list=None) -> xarray.Dataset:
    """
    reads FITS images and spatial az/el calibration for allsky camera
    Bdecl is in degrees, from IGRF model
    """
    warnings.filterwarnings('ignore', category=VerifyWarning)

    if isinstance(flist,str):
        flist = Path(flist).expanduser()

    if isinstance(flist,Path) and flist.is_dir():
        flist = list(flist.glob('*.FITS')) + list(flist.glob('*.fits'))
        flist = sorted(flist)

    flist = np.atleast_1d(flist)

    if treq is not None or wavelenreq is not None:
        time = []; wavelen=[]
        for i,fn in enumerate(flist):
            try:
                with fits.open(fn, mode='readonly') as h:
                    time.append(parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART']))
                    wavelen.append(int(h[0].header['FILTWAV']))
            except OSError: #many corrupted files, accounted for by preallocated vectors
                pass

        time = np.array(time)
        wavelen = np.array(wavelen)
#%% time request
    if treq is not None:
        if isinstance(treq,str):
            treq = parse(treq)
        elif isinstance(treq,(np.ndarray,tuple,list)) and isinstance(treq[0],str):
            treq = list(map(parse,treq))
        treq = np.atleast_1d(treq)
# %% time slice
        if treq.size == 1: # single frame
            i = abs(time-treq).argmin() #index number in flist desired
        elif treq.size == 2: #frames within bounds
            if treq[0] is not None:
                i = (treq[0] <= time)
                if treq[1] is not None:
                    i &= (time < treq[1])
            elif treq[1] is not None:
                i = (time < treq[1])
            else:
                i = slice(None)
        else:
            i = slice(None)

        if isinstance(i,np.ndarray) and not i.any():
            raise ValueError(f'no valid data found in {treq}')

        flist   = np.atleast_1d(flist[i])
        wavelen = np.atleast_1d(wavelen[i])
        if len(flist)==0:
            raise FileNotFoundError(f'no files found within time limits {treq}')
# %% wavelength slice
    if wavelenreq is not None:
        i = np.isin(wavelen, wavelenreq)
        flist = np.atleast_1d(flist[i])

    if len(flist)==0:
        raise FileNotFoundError(f'no files found with wavelength(s) {wavelenreq}')
#%% iterate over image files
    time = []; img= [];  wavelen = []
    for i,fn in enumerate(flist):
        try:
            with fits.open(fn, mode='readonly') as h:
                assert h[0].header['BITPIX']==16,'this function assumes unsigned 16-bit data'
                expstart = parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART'])

                time.append((expstart, expstart + timedelta(seconds=h[0].header['EXPTIME']))) #EXPTIME is in seconds

                wavelen.append(int(h[0].header['FILTWAV']))

                lla={'lat':h[0].header['GLAT'],
                     'lon':h[0].header['GLON'],
                     'alt_m':200.} # TODO use real altitude

                """
                DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
                just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
                since 14-bit \in {0,16384}.
                These files do not have a BZERO value. Someday when they're written correctly this
                code may need updating.
                Further, there was a RAID failure that filled the data files with random values.
                Don Hampton says about 90% of data OK, but 10% NOK.
                """

                I = np.rot90(h[0].data,-1) #NOTE: rotation to match online AVIs from UAF website. It's not transpose, and the cal file seems off.
                if not 'BZERO' in h[0].header.keys():
                    I[I>16384] = 0 #extreme, corrupted data
                    I = I.clip(0,16384).astype(np.uint16) #discard bad values for 14-bit cameras.

            img.append(I)

        except (OSError,TypeError) as e:
            logging.warning(f'{fn} has error {e}')

# %% collect output
    img = np.array(img)
    time = np.array(time)
    wavelen = np.array(wavelen)
    wavelengths = np.unique(wavelen)

    ds = {}
    for w in wavelengths:
        ds[w] = (('time','y','x'),img[wavelen==w,...])

    data = xarray.Dataset(ds,
                          coords={'time':time[:,0]},
                          attrs={'timeend':time[:,1],
                                 'lat':lla['lat'],'lon':lla['lon'],'alt_m':lla['alt_m']})

    if azfn is not None and elfn is not None:
        az,el = loadcal(azfn, elfn)
        data['az'] = az
        data['el'] = el

    return data


def loadcal(azfn:Path, elfn:Path) -> tuple:

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

    if (np.nanmax(el) > 90) or (np.nanmin(el) < 0):
        raise ValueError(' 0 < elevation < 90 degrees.')

    if (np.nanmax(az) > 360) or (np.nanmin(az) < 0):
        raise ValueError(' 0 < azimuth < 360 degrees.')

    el = (('y','x'),el)
    az = (('y','x'),az)

    return az,el