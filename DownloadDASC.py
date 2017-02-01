#!/usr/bin/env python
from pathlib import Path
from time import sleep
import ftplib
from pytz import UTC
from datetime import datetime
from dateutil.parser import parse
from urllib.parse import urlparse

def getdasc(start,end,host,site,odir='',clobber=False):
    """
    year,month,day: integer
    hour, minute:  start,stop integer len == 2
    """
    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2] + site
    odir = Path(odir).expanduser()
    print('downloading to {}'.format(odir.resolve()))
#%% get available files for this day
    rpath = fpath + '/DASC/RAW/{:4d}/{:4d}{:02d}{:02d}'.format(start.year,start.year,start.month,start.day)

    with ftplib.FTP(ftop,'anonymous','guest',timeout=15) as F:
        F.cwd(rpath)
        dlist = F.nlst()
        for f in dlist:
#%% file in time range
            #print (int(round(float(f[27:31]))))
            timg = datetime(int(f[14:18]), int(f[18:20]), int(f[20:22]),
                              int(f[23:25]), int(f[25:27]),tzinfo=UTC)
            if  start <= timg <= end:
#%% download file
                ofn = odir / f
                if not clobber:
                    if ofn.is_file(): #do filesizes match, if so, skip download
                        rsize = F.size(f)
                        if ofn.stat().st_size == rsize:
                            print('SKIPPING existing {}'.format(ofn))
                            continue

                print(ofn)
                with ofn.open('wb') as h:
                    F.retrbinary('RETR {}'.format(f), h.write)
                    sleep(1) # anti-leech


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('start',help='start time UTC e.g. 2012-11-03T06:23Z')
    p.add_argument('end',help='end time UTC e.g. 2012-11-03T06:25Z')
    p.add_argument('-o','--odir',help='directory to write downloaded FITS to',default='')
    p.add_argument('-c','--clobber',help='clobber (overwrite) existing files',action='store_true')
    p.add_argument('--host',default='ftp://optics.gi.alaska.edu')
    p.add_argument('-s','--site',help='EAA FYU KAK PKR TOO VEE',default='PKR')
    p = p.parse_args()

#host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"
    start = parse(p.start)
    end = parse(p.end)

    getdasc(start,end, p.host,p.site, p.odir, p.clobber)
