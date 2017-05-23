from dateutil.parser import parse
from datetime import datetime
from pytz import UTC
from numpy import ndarray,array

EPOCH = datetime(1970,1,1,tzinfo=UTC)

def totimestamp(t):
    """
    t may be None,float,int,datetime or list,tuple, ndarray of such
    output is ndarray of UTC Unix timestamp
    """

    if isinstance(t,datetime): # most cases devolve here
        t = t.timestamp()
    elif isinstance(t,str):
        t = totimestamp(parse(t))
    elif isinstance(t,(float,int)):
        t = float(t)
        assert 1e9 < t < 3e9, f'did you really mean {datetime.fromtimestamp(t,tz=UTC)}'
    elif isinstance(t,(tuple,list,ndarray)):
        t = array(map(totimestamp,t))

    return t

#print(totimestamp('2012-01-03T08:32:02Z'))
#print(totimestamp(1325579522))
#print(totimestamp(['2012-01-03T08:32:02Z','2012-01-03T08:32:12Z']))

def getdasc(start,end,host,site,odir='',clobber=False):
    """
    year,month,day: integer
    hour, minute:  start,stop integer len == 2
    """
    start = forceutc(start)
    end   = forceutc(end)

    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2] + site
    odir = Path(odir).expanduser()
    print('downloading to', odir.resolve())
#%% get available files for this day
    rpath = f'{fpath}/DASC/RAW/{start.year:4d}/{start.year:4d}{start.month:02d}{start.day:02d}'

    with ftplib.FTP(ftop,'anonymous','guest',timeout=15) as F:
        F.cwd(rpath)
        dlist = F.nlst()
        for f in dlist:
#%% file in time range
            #print (int(round(float(f[27:31]))))
            t = forceutc(datetime.strptime(f[14:-9],'%Y%m%d_%H%M%S'))
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
