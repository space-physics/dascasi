try:
    from pathlib import Path
    Path().expanduser()
except (ImportError,AttributeError):
    from pathlib2 import Path
#%%
from six import string_types,integer_types
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
        t = (t-EPOCH).total_seconds()
    elif isinstance(t,string_types):
        t = totimestamp(parse(t))
    elif isinstance(t,(float,integer_types)):
        t = float(t)
        assert 1e9 < t < 3e9, 'did you really mean {}'.format(datetime.fromtimestamp(t,tz=UTC))
    elif isinstance(t,(tuple,list,ndarray)):
        t = array([totimestamp(T) for T in t])

    return t

#print(totimestamp('2012-01-03T08:32:02Z'))
#print(totimestamp(1325579522))
#print(totimestamp(['2012-01-03T08:32:02Z','2012-01-03T08:32:12Z']))
