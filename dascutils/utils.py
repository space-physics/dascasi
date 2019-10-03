"""
Created on Fri Sep 28 11:43:24 2018

@author: smrak
"""
import xarray
import typing
from datetime import datetime, timedelta
from dateutil.parser import parse


def time_bounds(startend: typing.Tuple[datetime, datetime]) -> typing.Tuple[datetime, datetime]:
    start = parse(startend[0]) if isinstance(startend[0], str) else startend[0]  # type: ignore
    end = parse(startend[1]) if isinstance(startend[1], str) else startend[1]  # type: ignore

    if end < start:
        raise ValueError("start time must be before end time!")

    return start, end


def getPixelBrightness(
    D: xarray.Dataset, treq: typing.Tuple[datetime, datetime], obs_lat: float, obs_lon: float
) -> xarray.DataArray:
    """
    assuming fixed observation at obs_lat, obs_lon get brightness of a pixel_brightness

    for moving observer, loop this function.
    """
    # Filter according to time requirement
    i = (D.time >= treq[0]) & (D.time <= treq[-1])

    # GET BRIGHTNESS
    pixel_brightness = D.image[i].sel(lat=obs_lat, lon=obs_lon, method="nearest").values

    # Return brightness time-series
    return pixel_brightness


def getDASCimage(D: xarray.Dataset = None, ix: typing.Union[datetime, int] = None, coordinate: str = "wsg"):
    assert isinstance(ix, (datetime, int))
    # Find the closest image for the given timestamp
    if isinstance(ix, datetime):
        T = datetime2posix(ix)[0]
        Di = D.sel(time=T, method="nearest")
        dasc_dt = datetime.utcfromtimestamp(Di.time.values)
        if coordinate == "polar":
            img = Di.polar.values
        elif coordinate == "wsg":
            img = Di.image.values
    # Find for the given index
    if isinstance(ix, int):
        dasc_dt = datetime.utcfromtimestamp(D.time.values[ix])
        if coordinate == "polar":
            img = D.polar.values[ix]
        elif coordinate == "wsg":
            img = D.image.values[ix]

    return dasc_dt, img


def datetime2posix(dtime: typing.Union[datetime, typing.Sequence[datetime]]) -> typing.List[float]:
    """
    Convert an input list of datetime format timestamp to posix timestamp
    https://docs.python.org/3/library/datetime.html#datetime.datetime.timestamp
    """
    if isinstance(dtime, datetime):
        dtime = [dtime]

    return [(t - datetime(1970, 1, 1)) / timedelta(seconds=1) for t in dtime]
