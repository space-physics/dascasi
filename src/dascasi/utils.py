"""
Created on Fri Sep 28 11:43:24 2018

@author: smrak
"""

from datetime import datetime, timedelta
import numpy as np


def get_time_slice(time, treq: list[datetime] | None) -> slice:
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
        treq = [datetime.fromisoformat(treq)]
    if isinstance(treq[0], str):
        treq = [datetime.fromisoformat(treq[0]), datetime.fromisoformat(treq[1])]

    # %% time slice
    time = np.atleast_1d(time)  # type: ignore
    if len(treq) == 1:  # single frame
        j = abs(time - treq[0]).argmin()  # type: ignore
        i = slice(j, j + 1)  # ensures indexed lists remain list
    elif len(treq) == 2:  # frames within bounds
        i = slice(
            abs(time - treq[0]).argmin(),
            abs(time - treq[1]).argmin() + 1,  # type: ignore
        )  # type: ignore

    return i


def time_bounds(startend: tuple[datetime, datetime]) -> tuple[datetime, datetime]:
    start = (
        datetime.fromisoformat(startend[0])
        if isinstance(startend[0], str)
        else startend[0]
    )  # type: ignore
    end = (
        datetime.fromisoformat(startend[1])
        if isinstance(startend[1], str)
        else startend[1]
    )  # type: ignore

    if end < start:
        raise ValueError("start time must be before end time!")

    return start, end


def getDASCimage(D, ix: datetime | int, coordinate: str = "wsg"):
    # Find the closest image for the given timestamp
    match ix:
        case datetime():
            T = datetime2posix(ix)[0]
            Di = D.sel(time=T, method="nearest")
            dasc_dt = datetime.utcfromtimestamp(Di.time.values)
            if coordinate == "polar":
                img = Di.polar.values
            elif coordinate == "wsg":
                img = Di.image.values
        # Find for the given index
        case int():
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
