#!/usr/bin/env python
from dascutils import Path
import ftplib
from urllib.parse import urlparse

host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"

def getdasc(year,month,day,hour,minute,odir='',clobber=False):
    """
    year,month,day: integer
    hour, minute:  start,stop integer len == 2
    """
    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2]

#%% get available files for this day
    with ftplib.FTP(ftop,'anonymous','guest',timeout=15) as F:
        rpath = fpath + str(year) + '/{:4d}{:02d}{:02d}'.format(year,month,day)
        F.cwd(rpath)
        dlist = F.nlst()

        for f in dlist:
#%% file in time range?
            if ( hour[0] <= int(f[23:25]) <= hour[1] ) and ( minute[0] <= int(f[25:27]) <= minute[1] ):
#%% download file
                ofn = Path(odir).expanduser() / f
                if not clobber:
                    if ofn.is_file(): #do filesizes match, if so, skip download
                        rsize = F.size(f)
                        if ofn.stat().st_size == rsize:
                            print('SKIPPING locally existing {}'.format(ofn))
                            continue

                print(ofn)
                with ofn.open('wb') as h:
                    F.retrbinary('RETR {}'.format(f), h.write)


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
    P = p.parse_args()

    start = [int(i) for i in P.start.split(':')]
    end =   [int(i) for i in P.end.split(':')]

    hour = (start[0],end[0])
    minute = (start[1],end[1])

    getdasc(P.year, P.month, P.day, hour, minute, P.odir,P.clobber)