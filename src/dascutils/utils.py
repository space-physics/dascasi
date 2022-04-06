"""
Created on Fri Sep 28 11:43:24 2018

@author: smrak
"""

from __future__ import annotations
import xarray
from datetime import datetime, timedelta
from dateutil.parser import parse
import numpy as np


def get_time_slice(time, treq: list[datetime]) -> slice:
    """
    given the times in a data stack and the requested time(s),
    return a slice for indexing the data stack

    Parameters
    ----------
    time : list of datetime
        times in data stack
    treq : list of datetime
        requested time(s)

    Returns
    -------
    i: slice
        indices corresponding to requested time(s) in data stack

    """
    if treq is None:
        return slice(None)
    if isinstance(treq, str):
        treq = [parse(treq)]
    if isinstance(treq[0], str):
        treq = [parse(treq[0]), parse(treq[1])]

    # %% time slice
    time = np.atleast_1d(time)  # type: ignore
    if len(treq) == 1:  # single frame
        j = abs(time - treq[0]).argmin()  # type: ignore
        i = slice(j, j + 1)  # ensures indexed lists remain list
    elif len(treq) == 2:  # frames within bounds
        i = slice(
            abs(time - treq[0]).argmin(), abs(time - treq[1]).argmin() + 1  # type: ignore
        )  # type: ignore

    return i


def time_bounds(startend: tuple[datetime, datetime]) -> tuple[datetime, datetime]:
    start = parse(startend[0]) if isinstance(startend[0], str) else startend[0]  # type: ignore
    end = parse(startend[1]) if isinstance(startend[1], str) else startend[1]  # type: ignore

    if end < start:
        raise ValueError("start time must be before end time!")

    return start, end


def getDASCimage(D: xarray.DataArray, ix: datetime | int = None, coordinate: str = "wsg"):
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


def datetime2posix(dtime: datetime | list[datetime]) -> list[float]:
    """
    Convert an input list of datetime format timestamp to posix timestamp
    https://docs.python.org/3/library/datetime.html#datetime.datetime.timestamp
    """
    if isinstance(dtime, datetime):
        dtime = [dtime]

    return [(t - datetime(1970, 1, 1)) / timedelta(seconds=1) for t in dtime]
