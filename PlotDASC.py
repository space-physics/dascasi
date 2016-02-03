#!/usr/bin/env python3
"""
Plots / plays Poker Flat DASC all-sky camera data FITS files

To download DASC images using Octave, Matlab, or Python checkout:
https://github.com/jswoboda/ISR_Toolbox/blob/master/Example_Scripts/loadDASC2013Apr14.py
or
download manually from
https://amisr.asf.alaska.edu/PKR/DASC/RAW/
note the capitalization is required in that URL.
"""
from dascutils.readDASCfits import readallDasc
from dascutils.plotdasc import histdasc,moviedasc
#
from themisasi.readthemis import mergefov

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='for Poker Flat DASC all sky camera, read az/el mapping and images')
    p.add_argument('indir',help='directory of .fits or specific .fits file')
    p.add_argument('-a','--azfn',help='filename for DASC .fits azimuth calibration')
    p.add_argument('-e','--elfn',help='filename for DASC .fits elevation calibration')
    p.add_argument('-w','--wavelength',help='select wavelength(s) to plot simultaneously [428 558 630]',type=int,default=[428,558,630],nargs='+')
    p.add_argument('-m','--minmax',help='set values outside these limits to 0, due to data corruption',type=int,nargs=2,default=[350,9000])
    p.add_argument('-c','--cadence',help='set playback cadence to request times [sec]',type=float,default=5.)
    p.add_argument('-o','--odir',help='output directory')
    p.add_argument('--ncal',help='narrow FOV camera calibration files HDF5',nargs='+')
    p.add_argument('--projalt',help='altitude [km] to project common FOV at',type=float,default=110)
    p=p.parse_args()

    ocalfn = None  #filename to save az,el contour plot png

    img,times,waz,wel,wlla = readallDasc(p.indir,p.azfn,p.elfn,p.wavelength,p.minmax)

    rows,cols = mergefov(ocalfn,wlla,waz,wel,None,None,p.ncal,p.projalt)

#%% plots
    histdasc(img,p.wavelength,p.odir)

    moviedasc(img,p.wavelength,times,p.odir,p.cadence)
