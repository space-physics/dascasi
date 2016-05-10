#!/usr/bin/env python3
"""
Reads DASC allsky cameras images in FITS formats into GeoData.
Run standalone from PlayDASC.py
"""
import logging
from pathlib import Path
from astropy.io import fits
import numpy as np
from dateutil.parser import parse
from warnings import warn
#
from histutils.fortrandates import forceutc

def readallDasc(indir,azfn,elfn,wl,minmax):
    """
    returns Dasc images in list by wavelength, then by time
    """

    img = []; times = []
    for w in wl:
        data,coordnames,dataloc,sensorloc,time  = readCalFITS(indir,azfn,elfn,w,minmax)
        img.append(data['image'])
        times.append(time)
#%% histogram
    try:
        az = dataloc[:,1].reshape(img[0].shape[1:])
        el = dataloc[:,2].reshape(img[0].shape[1:])
    except Exception: #azel data wasn't loaded
        az=el=None

    return img,times,az,el,sensorloc

def readCalFITS(indir,azfn,elfn,wl,minmax):
    indir = Path(indir).expanduser()
    if not wl:
        wl = '*' #select all wavelengths

    flist = []
    #for w in wl:
    flist += sorted(indir.glob("PKR_DASC_0{}_*.FITS".format(wl)))
    return readDASC(flist,azfn,elfn,minmax=minmax)

def readDASC(flist,azfn,elfn,heightkm=110,minmax=None):
    """
    reads FITS images and spatial az/el calibration for allsky camera
    """
    if not flist:
        warn('no data files found')
        return
#%% preallocate, assuming all images the same size
    with fits.open(str(flist[0]),mode='readonly') as h:
        img = h[0].data
    sensorloc = np.empty(3) #in case error in reading this file
    dataloc = np.empty((img.size,3))
    times =   np.empty((len(flist),2)); times.fill(np.nan)
    assert h[0].header['BITPIX']==16,'this function assumes unsigned 16-bit data'
    img =     np.zeros((len(flist),img.shape[0],img.shape[1]),np.uint16) #zeros in case a few images fail to load
    wavelen = np.empty(len(flist)); wavelen.fill(np.nan)
#%% iterate over image files
    for i,fn in enumerate(flist):
        try:
            with fits.open(str(fn),mode='readonly') as h:
                expstart_dt = forceutc(parse(h[0].header['OBSDATE'] + ' ' + h[0].header['OBSSTART']))
                expstart_unix = expstart_dt.timestamp()
                times[i,:] = [expstart_unix,expstart_unix + h[0].header['EXPTIME']]

                wavelen[i] = h[0].header['FILTWAV']

                sensorloc={'lat':h[0].header['GLAT'],
                           'lon':h[0].header['GLON'],
                            'alt_m':200.} #TODO use real DASC altitude

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
                    I = I.clip(0,16384).astype(np.uint16) #discard bad values for 14-bit cameras.
                    I = np.rot90(I,-1) #NOTE: rotation to match online AVIs from UAF website

                img[i,...] = I

        except Exception as e:
            logging.info('{} has error {}'.format(fn,e))


    if minmax is not None:
        #%% deal with corrupted data
        img[(img<minmax[0]) | (img>minmax[1])] = 1 #instead of 0 for lognorm

    #we return the images as a 3-D array data['image'] all wavelengths stacked together, which you can select by data['lambda']
    data = {'image':img,
            'lambda':wavelen}

    coordnames="spherical"
    try:
        with fits.open(str(Path(azfn).expanduser()),mode='readonly') as h:
            az = np.rot90(h[0].data,2) # NOTE: rotation to match UAF AVIs NOT flipud
        with fits.open(str(Path(elfn).expanduser()),mode='readonly') as h:
            el = np.rot90(h[0].data,2) # NOTE: rotation to match UAF AVIs NOT flipud
        dataloc[:,0] = heightkm
        dataloc[:,1] = az.ravel()
        dataloc[:,2] = el.ravel()
    except Exception as e:
       dataloc=None

    return data,coordnames,dataloc,sensorloc,times