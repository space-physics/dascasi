#!/usr/bin/env python
"""
Plots / plays / converts to movie:  Poker Flat DASC all-sky camera data FITS files

example:
python PlotDASC.py tests/PKR_DASC_0428_20151007_082305.930.FITS
azfn = R.parent/'cal/PKR_DASC_20110112_AZ_10deg.fits
elfn = R.parent/'cal/PKR_DASC_20110112_EL_10deg.fits
"""
from pathlib import Path
import xarray
from matplotlib.pyplot import show
from argparse import ArgumentParser
import dascutils.io as dio
import dascutils.plots as dup


def plotdasc(img: xarray.Dataset, odir: Path, cadence: float):
    dup.histogram_dasc(img, odir)

    dup.moviedasc(img, odir, cadence)


def main():
    p = ArgumentParser(description='for Poker Flat DASC all sky camera, read az/el mapping and images')
    p.add_argument('indir', help='directory of .fits or specific .fits file')
    p.add_argument('-t', '--tlim', help='only plot data in this range', nargs=2)
    p.add_argument('-a', '--azelfn', help='stem for DASC .fits azimuth calibration',
                   default='cal/PKR_DASC_20110112')
    p.add_argument('-w', '--wavelength', help='select wavelength(s) to plot simultaneously [428 558 630]', type=int, nargs='+')
    p.add_argument('-m', '--minmax', help='set values outside these limits to 0, due to data corruption',
                   type=int, nargs=2, default=[350, 9000])
    p.add_argument('-c', '--cadence', help='set playback cadence to request times [sec]', type=float, default=5.)
    p.add_argument('-o', '--odir', help='output directory')
    p = p.parse_args()

    imgs = dio.load(p.indir, p.azelfn, p.tlim, p.wavelength)

    plotdasc(imgs, p.odir, p.cadence)

    show()


if __name__ == '__main__':
    main()
