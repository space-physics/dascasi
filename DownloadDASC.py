#!/usr/bin/env python
from dateutil.parser import parse
#
from dascutils import getdasc

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('start',help='start time UTC e.g. 2012-11-03T06:23Z')
    p.add_argument('end',help='end time UTC e.g. 2012-11-03T06:25Z')
    p.add_argument('-o','--odir',help='directory to write downloaded FITS to',default='')
    p.add_argument('-c','--clobber',help='clobber (overwrite) existing files',action='store_true')
    p.add_argument('-host',default='ftp://optics.gi.alaska.edu')
    p.add_argument('-s','--site',help='EAA FYU KAK PKR TOO VEE',default='PKR')
    p = p.parse_args()

#host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"
    start = parse(p.start)
    end = parse(p.end)

    getdasc(start,end, p.host,p.site, p.odir, p.clobber)
