from pathlib import Path
from time import sleep
import ftplib
from dateutil.parser import parse
from datetime import datetime
from urllib.parse import urlparse


def totimestamp(t):
    """
    t may be None,float,int,datetime or list,tuple, ndarray of such
    output is ndarray of UTC Unix timestamp
    """
    if t is None:
        return

    if isinstance(t,datetime): # most cases devolve here
        t = t.timestamp()
    elif isinstance(t,str):
        t = totimestamp(parse(t))
    elif isinstance(t,(float,int)):
        t = float(t)
        assert 1e9 < t < 3e9, f'did you  mean {datetime.fromtimestamp(t)}'
    else: # assume it's an iterable 1-D vector
        t = list(map(totimestamp,t))

    return t


def download(startend,host,site,odir='',clobber=False):
    """
    startend: tuple of datetime
    year,month,day: integer
    hour, minute:  start,stop integer len == 2
    """
    assert len(startend)==2

    start = parse(startend[0]) if isinstance(startend[0],str) else startend[0]
    end   = parse(startend[1]) if isinstance(startend[1],str) else startend[1]

    assert end >= start,'start time must be before end time!'

    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2] + site
    odir = Path(odir).expanduser().resolve()
    if not odir.is_dir():
        raise FileNotFoundError(f'{odir} does not exist')
#%% get available files for this day
    rparent = f'{fpath}/DASC/RAW/{start.year:4d}'
    rday = f'{start.year:4d}{start.month:02d}{start.day:02d}'

    with ftplib.FTP(ftop,'anonymous','guest',timeout=15) as F:
        F.cwd(rparent)
        dlist = F.nlst()
        if not rday in dlist:
            raise FileNotFoundError(f'{rday} does not exist under {host}/{rparent}')

        print('downloading to', odir)
        F.cwd(rday)
        dlist = F.nlst()

        print(f'remote filesize approx. {F.size(dlist[0])/1000} kB.')
        for f in dlist:
#%% file in time range
            #print (int(round(float(f[27:31]))))
            t = datetime.strptime(f[14:-9],'%Y%m%d_%H%M%S')
            if  start <= t <= end:
#%% download file
                ofn = odir / f
                if not clobber:
                    if ofn.is_file(): #do filesizes match, if so, skip download
                        rsize = F.size(f)
                        if ofn.stat().st_size == rsize:
                            print('SKIPPING existing', ofn)
                            continue

                print(ofn)
                with ofn.open('wb') as h:
                    F.retrbinary(f'RETR {f}', h.write)
                    sleep(1) # anti-leech
