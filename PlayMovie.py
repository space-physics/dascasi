#!/usr/bin/env python
"""
Plots / plays / converts to movie:  Poker Flat DASC all-sky camera data FITS files

example:
python PlayMovie.py tests/
"""
from pathlib import Path
import xarray
from matplotlib.pyplot import show
from argparse import ArgumentParser
import dascutils as du
import dascutils.plots as dup

R = Path(__file__).parent.resolve()


def plotdasc(img: xarray.Dataset, outdir: Path, cadence: float):
    dup.histogram_dasc(img, outdir)

    dup.moviedasc(img, outdir, cadence)


def main():
    p = ArgumentParser(description="for DASC all sky camera, read and play az/el mapping and images")
    p.add_argument("indir", help="directory of .fits or specific .fits file")
    p.add_argument("-o", "--outdir", help="directory to write plots to")
    p.add_argument("-t", "--tlim", help="only plot data in this range", nargs=2)
    p.add_argument("-a", "--azelfn", help="stem for DASC .fits azimuth calibration", default=R / "cal/PKR_DASC_20110112")
    p.add_argument("-w", "--wavelength", help="select wavelength(s) to plot simultaneously [428 558 630]", nargs="+")
    p.add_argument("-c", "--cadence", help="set playback cadence to request times [sec]", type=float, default=5.0)
    p = p.parse_args()

    imgs = du.load(p.indir, p.azelfn, treq=p.tlim, wavelenreq=p.wavelength)

    plotdasc(imgs, p.outdir, p.cadence)

    show()


if __name__ == "__main__":
    main()
