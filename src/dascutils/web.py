#!/usr/bin/env python
"""
I tried to parallelize this via threading etc, but the FTP server just trips the anti-leech
and sends a bunch of zero-sized files.

Best to download once and convert FITS stack to HDF5.
"""

from __future__ import annotations
import typing as T
from pathlib import Path
from time import sleep
import ftplib
from datetime import datetime
from urllib.parse import urlparse

from .utils import time_bounds

HOST = "ftp://optics.gi.alaska.edu"


def download(
    startend: tuple[datetime, datetime],
    site: str,
    odir: Path,
    host: str = None,
    wavelen: str = None,
) -> list[Path]:
    """
    startend: tuple of datetime
    """
    if not host:
        host = HOST

    assert len(startend) == 2

    start, end = time_bounds(startend)

    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2] + site
    odir = Path(odir).expanduser().resolve()
    odir.mkdir(exist_ok=True, parents=True)
    # %% get available files for this day
    rparent = f"{fpath}/DASC/RAW/"
    rday = f"{start.year:4d}{start.month:02d}{start.day:02d}"
    # %% wavelength
    if wavelen is None:
        pass
    elif isinstance(wavelen, int):
        wavelen = f"{wavelen:04d}"
    elif isinstance(wavelen, str):
        if len(wavelen) != 4:
            raise ValueError("expecting 4-character wavelength spec e.g. 0428")
    elif not isinstance(wavelen, (tuple, list)):
        raise TypeError("expecting 4-character wavelength spec e.g. 0428")

    flist = []

    with ftplib.FTP(ftop, "anonymous", "guest", timeout=15) as F:
        slist = F.nlst()
        if site not in slist:
            raise ValueError(f"site '{site}' not found in remote directory: {slist}")
        F.cwd(rparent)
        ylist = F.nlst()
        if str(start.year) not in ylist:
            raise ValueError(f"site {site} only has data for years: {ylist}")
        ydir = str(start.year)
        F.cwd(ydir)
        dlist = F.nlst()
        if rday not in dlist:
            raise ValueError(
                f"{rday} does not exist under {host}/{rparent}/{ydir}. Available: {dlist}"
            )

        print("downloading to", odir)
        F.cwd(rday)

        for filename in get_filenames(F.nlst(), wavelen, start, end):
            # %% download file
            ofn = odir / filename
            flist.append(ofn)

            if skip_exist(ofn, F):
                continue

            print(ofn)
            with ofn.open("wb") as h:
                F.retrbinary(f"RETR {filename}", h.write)
                sleep(0.5)  # anti-leech

        if len(flist) == 0:
            raise ValueError(f"could not find files between {start} {end}")

        print(f"found or downloaded {len(flist)} files")

    return flist


def skip_exist(filename: Path, F) -> bool:
    if not filename.is_file():
        return False
    if filename.stat().st_size == F.size(filename.name):
        print("SKIPPING existing", filename)
        return True
    return False


def get_filenames(days: list[str], wavelen: str, start: datetime, end: datetime) -> T.Iterator[str]:

    print(f"searching {len(days)} remote files")

    for filename in days:
        # %% qualifiers
        if not filename.lower().endswith((".fit", ".fits")):
            # extra subdirs, extraneous files
            continue
        if wavelen and filename[9:13] not in wavelen:
            continue

        # names like:
        # * VEE_DASC_0000_20170102_030405.FIT
        # * PKR_DASC_0428_20170102_030405.678.FITS
        tfile = datetime.strptime(filename[14:29], "%Y%m%d_%H%M%S")
        if tfile < start or tfile > end:
            continue
        yield filename
