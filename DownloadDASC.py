#!/usr/bin/env python
import dascutils as du
from argparse import ArgumentParser


def main():
    p = ArgumentParser(description="download DASC all-sky camera data")
    p.add_argument("site", help="EAA FYU KAK PKR TOO VEE")
    p.add_argument("startend", help="start/end times UTC e.g. 2012-11-03T06:23 2012-11-03T07", nargs=2)
    p.add_argument("odir", help="directory to write downloaded FITS to")
    p.add_argument("-w", "--wavelen", help="request specific wavelength(s)", nargs="+")
    p.add_argument("-host", default="ftp://optics.gi.alaska.edu")
    p = p.parse_args()

    # host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"
    du.download(p.startend, p.site, p.odir, p.host, p.wavelen)


if __name__ == "__main__":
    main()
