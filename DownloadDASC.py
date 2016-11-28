#!/usr/bin/env python
from dascutils import Path
from time import sleep
import ftplib
from datetime import datetime
from six.moves.urllib.parse import urlparse

def getdasc(year,month,day,hour,minute,host,site,odir='',clobber=False):
    """
    year,month,day: integer
    hour, minute:  start,stop integer len == 2
    """
    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2] + site

#%% get available files for this day
    with ftplib.FTP(ftop,'anonymous','guest',timeout=15) as F:
        rpath = fpath + '/DASC/RAW/{:4d}/{:4d}{:02d}{:02d}'.format(year,year,month,day)
        F.cwd(rpath)
        dlist = F.nlst()
        start = datetime(year, month, day, hour[0], minute[0]).timestamp()
        stop = datetime(year, month, day, hour[1], minute[1]).timestamp()
        for f in dlist:
#%% file in time range
            #print (int(round(float(f[27:31]))))
            img_ts = datetime(int(f[14:18]), int(f[18:20]), int(f[20:22]),
                              int(f[23:25]), int(f[25:27]), int(round(float(f[27:31])))%60).timestamp()
            if (img_ts >= start) and (img_ts <= stop):
#%% download file
                ofn = Path(odir).expanduser() / f
                print(ofn)
                if not clobber:
                    if ofn.is_file(): #do filesizes match, if so, skip download
                        rsize = F.size(f)
                        if ofn.stat().st_size == rsize:
                            print('SKIPPING locally existing {}'.format(ofn))
                            continue

                print(ofn)
                with ofn.open('wb') as h:
                    F.retrbinary('RETR {}'.format(f), h.write)
                    sleep(1)


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('year',type=int)
    p.add_argument('month',type=int)
    p.add_argument('day',type=int)
    p.add_argument('start',help='start time UTC e.g. 06:23')
    p.add_argument('end',help='end time UTC e.g. 06:25')
    p.add_argument('-o','--odir',help='directory to write downloaded FITS to',default='')
    p.add_argument('-c','--clobber',help='clobber (overwrite) existing files',action='store_true')
    p.add_argument('--host',default='ftp://optics.gi.alaska.edu')
    p.add_argument('-s','--site',help='EAA FYU KAK PKR TOO VEE',default='PKR')
    P = p.parse_args()

#host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"

    start = [int(i) for i in P.start.split(':')]
    end =   [int(i) for i in P.end.split(':')]

    hour = (start[0],end[0])
    minute = (start[1],end[1])

    getdasc(P.year, P.month, P.day, hour, minute, P.host,P.site, P.odir,P.clobber)
