#!/usr/bin/env python
"""
Reads DASC allsky cameras images in FITS formats into GeoData.
Run standalone from PlayDASC.py
"""
from pathlib import Path
from warnings import filterwarnings # corrupt FITS files let off a flood of AstroPy warnings
from astropy.io.fits.verify import VerifyWarning
import logging
from astropy.io import fits
import numpy as np
from dateutil.parser import parse
#
from . import totimestamp
#
from sciencedates import forceutc


def readallDasc(indir,azfn,elfn,wl,minmax,tlim=None):
    """
    returns Dasc images in list by wavelength, then by time
    """

    img = []; times = []; wlused=[]
    for w in wl:
        try:
            data,azel,sensorloc,time  = readCalFITS(indir,azfn,elfn,w,minmax,tlim)
            img.append(data['image'])
            times.append(time)
            wlused.append(w)
        except FileNotFoundError:
            pass
#%% histogram
    try:
        az = azel[0]
        el = azel[1]
    except Exception: #azel data wasn't loaded
        az=el=None

    return img,times,az,el,sensorloc,wlused

def readCalFITS(indir,azfn,elfn,wl,minmax,tlim=None):
    indir = Path(indir).expanduser()
    if not wl:
        wl = '*' #select all wavelengths

    #flist = []
    #for w in wl:
    flist = sorted(indir.glob("PKR_DASC_0{}_*.FITS".format(wl)))
    return readDASC(flist,azfn,elfn,minmax,tlim)

def readDASC(flist,azfn=None,elfn=None,minmax=None,treq=None):
    """
    reads FITS images and spatial az/el calibration for allsky camera
    Bdecl is in degrees, from IGRF model
    """
    filterwarnings('ignore',category=VerifyWarning)

    if not flist:
        raise FileNotFoundError('no files of this wavelength')

    flist = np.atleast_1d(flist)

    treq = totimestamp(treq)
#%% read one file mode
    if treq is not None:
        expstart = np.empty(len(flist)); expstart.fill(np.nan)

        for i,fn in enumerate(flist):
            try:
                with fits.open(str(fn),mode='readonly') as h:
                    expstart[i] = forceutc(parse(h[0].header['OBSDATE'] + ' ' + h[0].header['OBSSTART'])).timestamp()
            except IOError: #many corrupted files, accounted for by preallocated vectors
                pass



        if isinstance(treq,float) or len(treq) == 1: # single frame
            fi = np.nanargmin(abs(expstart-treq)) #index number in flist desired
        elif len(treq)==2: #frames within bounds
            if treq[0] is not None:
                fi = (treq[0] <= expstart)
                if treq[1] is not None:
                    fi &= (expstart < treq[1])
            elif treq[1] is not None:
                fi = (expstart < treq[1])
            else:
                fi = slice(None)
        else:
            fi = slice(None)

        flist = flist[fi]
        if len(flist)==0:
            raise FileNotFoundError('no files found within time limits')

        if isinstance(flist,Path): # so that we can iterate
            flist = [flist]

#%% preallocate, assuming all images the same size
    for f in flist: #find the first "good" file
        try:
            with fits.open(str(f),mode='readonly') as h:
                img = h[0].data
                sensorloc={'lat':h[0].header['GLAT'],
                           'lon':h[0].header['GLON'],
                            'alt_m':200.} #TODO use real DASC altitude
            break
        except IOError: # taking first good file
            pass

    times =   np.empty((len(flist),2)); times.fill(np.nan)
    assert h[0].header['BITPIX']==16,'this function assumes unsigned 16-bit data'
    img =     np.zeros((len(flist),img.shape[0],img.shape[1]),np.uint16) #zeros in case a few images fail to load
    wavelen = np.empty(len(flist)); wavelen.fill(np.nan)
    iok = np.zeros(len(flist),bool)
#%% iterate over image files
    for i,fn in enumerate(flist):
        try:
            with fits.open(str(fn),mode='readonly') as h:
                expstart = forceutc(parse(h[0].header['OBSDATE'] + ' ' + h[0].header['OBSSTART'])).timestamp()

                times[i,:] = [expstart,expstart + h[0].header['EXPTIME']] #EXPTIME is in seconds

                wavelen[i] = h[0].header['FILTWAV']

                """
                DASC iKon cameras are/were 14-bit at least through 2015. So what they did was
                just write unsigned 14-bit data into signed 16-bit integers, which doesn't overflow
                since 14-bit \in {0,16384}.
                These files do not have a BZERO value. Someday when they're written correctly this
                code may need updating.
                Further, there was a RAID failure that filled the data files with random values.
                Don Hampton says about 90% of data OK, but 10% NOK.
                """

                I = h[0].data
                if not 'BZERO' in h[0].header.keys():
                    I[I>16384] = 0 #extreme, corrupted data
                    I = I.clip(0,16384).astype(np.uint16) #discard bad values for 14-bit cameras.

            img[i,...] =  np.rot90(I,-1) #NOTE: rotation to match online AVIs from UAF website. It's not transpose, and the cal file seems off.
            iok[i] = True

        except (IOError,TypeError) as e:
            logging.info('{} has error {}'.format(fn,e))

#%% keep only good times
    img = img[iok,...]
    times = times[iok,:]
    wavelen = wavelen[iok]
#%% deal with corrupted data
    if minmax is not None:
        img[(img<minmax[0]) | (img>minmax[1])] = 1 #instead of 0 for lognorm

    #we return the images as a 3-D array data['image'] all wavelengths stacked together, which you can select by data['lambda']
    data = {'image':img,
            'lambda':wavelen}

    if azfn is not None and elfn is not None:
        with fits.open(str(Path(azfn).expanduser()),mode='readonly') as h:
            az = h[0].data
        with fits.open(str(Path(elfn).expanduser()),mode='readonly') as h:
            el = h[0].data # NOTE: no rotation/flip
    else:
        az=el=None

    return data,(az,el),sensorloc,times
