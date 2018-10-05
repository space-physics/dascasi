#!/usr/bin/env python
import dascutils as du
from argparse import ArgumentParser


def main():
    p = ArgumentParser()
    p.add_argument('startend', help='start/end times UTC e.g. 2012-11-03T06:23', nargs=2)
    p.add_argument('odir', help='directory to write downloaded FITS to')
    p.add_argument('-c', '--overwrite', help='overwrite existing files', action='store_true')
    p.add_argument('-host', default='ftp://optics.gi.alaska.edu')
    p.add_argument('-s', '--site', help='EAA FYU KAK PKR TOO VEE', default='PKR')
    p.add_argument('-w', '--wavelen', help='request specific wavelength(s)', nargs='+')
    p = p.parse_args()

    # host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"
    du.download(p.startend, p.site, p.odir, p.host, p.wavelen, p.overwrite)


if __name__ == '__main__':
    main()
