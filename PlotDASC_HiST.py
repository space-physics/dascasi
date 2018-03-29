#!/usr/bin/env python
"""
Plots / plays / converts to movie:  Poker Flat DASC all-sky camera data FITS files

This program by default projects HiST auroral tomography system FOV onto PFRR DASC.
"""
import matplotlib
#matplotlib.use('Agg')
#
from dascutils.readDASCfits import readallDasc
from dascutils.plots import histdasc,moviedasc
#
from themisasi.fov import mergefov

def plothstfovondasc(img,wavelength,odir,cadence,rows,cols):
    histdasc(img,wavelength,odir) #histogram

    moviedasc(img,wavelength,times,odir,cadence,rows,cols)


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='for Poker Flat DASC all sky camera, read az/el mapping and images')
    p.add_argument('indir',help='directory of .fits or specific .fits file')
    p.add_argument('-t','--tlim',help='only plot data in this range',nargs=2,default=(None,None))
    p.add_argument('-a','--azfn',help='filename for DASC .fits azimuth calibration',default='cal/PKR_DASC_20110112_AZ_10deg.fits')
    p.add_argument('-e','--elfn',help='filename for DASC .fits elevation calibration',default='cal/PKR_DASC_20110112_EL_10deg.fits')
    p.add_argument('-w','--wavelength',help='select wavelength(s) to plot simultaneously [428 558 630]',type=int,default=[428,558,630],nargs='+')
    p.add_argument('-m','--minmax',help='set values outside these limits to 0, due to data corruption',type=int,nargs=2,default=[350,9000])
    p.add_argument('-c','--cadence',help='set playback cadence to request times [sec]',type=float,default=5.)
    p.add_argument('-o','--odir',help='output directory',default='.')
    p.add_argument('--ncal',help='narrow FOV camera calibration files HDF5',nargs='+',default=['../histutils/cal/hst0cal.h5','../histutils/cal/hst1cal.h5'])
    p.add_argument('--projalt',help='altitude [METERS] to project common FOV at',type=float,default=110e3)
    p=p.parse_args()

    ocalfn = None  #filename to save az,el contour plot png

    try:
        plothstfovondasc(img,p.wavelength,p.odir,p.cadence,rows,cols)
    except NameError:
        img,times,waz,wel,wlla = readallDasc(p.indir,p.azfn,p.elfn,p.wavelength,p.minmax,p.tlim)
        rows,cols = mergefov(ocalfn,wlla,waz,wel,None,None,p.ncal,p.projalt,site='DASC')

        plothstfovondasc(img,p.wavelength,p.odir,p.cadence,rows,cols)
