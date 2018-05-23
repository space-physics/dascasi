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
         wavelenreq:list=None, verbose=False) -> xarray.Dataset:
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

    flist = np.unique(flist)
    if treq is not None or wavelenreq is not None:
        time = []; wavelen=[]; flist1=[]
        for i,fn in enumerate(flist):
            try:
                with fits.open(fn, mode='readonly') as h:
                    # New DASC
                    if 'OBSDATE' in h[0].header:
                        time.append(parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART']))
                        wavelen.append(int(h[0].header['FILTWAV']))
                        flist1.append(fn)
                    # Old DASC has differnent header
                    else:
                        time.append(parse(h[0].header['FRAME']))
                        wavelen.append(0)
                        flist1.append(fn)
            except OSError: #many corrupted files, accounted for by preallocated vectors
                pass
        flist = np.array(flist1)
        time = np.array(time)
        wavelen = np.array(wavelen)
        print ('Wavelenghts detected in the folder: {}'.format(np.unique(wavelen)))
#%% time request
    if treq is not None:
        if isinstance(treq,str):
            treq = parse(treq)
        elif isinstance(treq,(np.ndarray,tuple,list)) and isinstance(treq[0],str):
            treq = list(map(parse,treq))
        elif isinstance(treq,(np.ndarray,tuple,list)) and isinstance(treq[0],datetime):
            pass
        if isinstance(treq,list):
            treq = np.array(treq)
# %% time slice
        if treq.size == 1: # single frame
            i = abs(time-treq).argmin() #index number in flist desired
        elif treq.size == 2: #frames within bounds
            if treq[0] is not None:
                i = np.where( (time>=treq[0]) & (time<=treq[1]) )[0]

        if isinstance(i,np.ndarray) and not i.any():
            raise ValueError(f'no valid data found in {treq}')
        flist = flist[i]
        wavelen = wavelen[i]
        if len(flist)==0:
            raise FileNotFoundError(f'no files found within time limits {treq}')
# %% wavelength slice
    if wavelenreq is not None:
        i = np.isin(wavelen, wavelenreq)
        flist = flist[i]

    if len(flist)==0:
        raise FileNotFoundError(f'no files found with wavelength(s) {wavelenreq}')
#%% iterate over image files
    time = []; img= [];  wavelen = []
    for i,fn in enumerate(flist):
        try:
            with fits.open(fn, mode='readonly', memmap=False) as h:
                if verbose and i == 0:
                    print ('Number of files: ',len(flist))
                    print (h[0].header)
                assert h[0].header['BITPIX']==16,'this function assumes unsigned 16-bit data'
                if 'OBSDATE' in h[0].header and 'EXPTIME' in h[0].header:
                     time.append(parse(h[0].header['OBSDATE'] + 'T' + h[0].header['OBSSTART']))
                elif 'FRAME' in h[0].header: #old DASC files
                    time.append(parse(h[0].header['FRAME']))
                    wavelen.append(0)
                try:
                    wavelen.append(int(h[0].header['FILTWAV']))
                except KeyError:
                    pass

                try:
                    lla={'lat':h[0].header['GLAT'],
                         'lon':h[0].header['GLON'],
                         'alt_m':200.} # TODO use real altitude
                except KeyError:
                    # Hard coded for Poker Flat
                    lla = {'lat':65.12992,
                         'lon':-147.47104,
                         'alt_m':200.}
                """
                DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
                just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
                since 14-bit \in {0,16384}.
                These files do not have a BZERO value. Someday when they're written correctly this
                code may need updating.
                Further, there was a RAID failure that filled the data files with random values.
                Don Hampton says about 90% of data OK, but 10% NOK.
                """
                # Az == 0 toward east. Rot90 to match the north=0deg.
                I = np.rot90(h[0].data) #NOTE: to match the online AVIs from UAF website, rotate the file with np.rot90(I,-1).
                if not 'BZERO' in h[0].header:
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
    if wavelengths.shape[0] == 1 and wavelengths[0] == 0:
        img = img[:,:,0,:]
        ds[0] = (('time','y','x'),img[wavelen==0,...])
    else:
        for w in wavelengths:
            print (img[wavelen==w,...].shape)
            ds[w] = (('time','y','x'),img[wavelen==w,...])

    data = xarray.Dataset(ds,
                          coords={'time':time},)
    if lla is not None:
        data.attrs['lat']=lla['lat']
        data.attrs['lon']=lla['lon']
        data.attrs['alt_m']=lla['alt_m']

    if azfn is not None and elfn is not None:
        az,el = loadcal(azfn, elfn)
        data.attrs['az']=az
        data.attrs['el']=el

    data.attrs['filename'] = ' '.join((p.name for p in flist))

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
    bad &= el==0

    el[bad] = np.nan
    az[bad] = np.nan

    assert np.nanmax(el) <= 90  and np.nanmin(el) >= 0, '0 < elevation < 90 degrees.'
    assert np.nanmax(az) <= 360 and np.nanmin(az) >= 0, '0 < azimuth < 360 degrees.'

    el = (('y','x'), el)
    az = (('y','x'), az)

    return az,el