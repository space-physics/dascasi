#!/usr/bin/env python
from pathlib import Path
from time import sleep
import ftplib
from dateutil.parser import parse
from datetime import datetime
from typing import Tuple, List
from urllib.parse import urlparse

HOST = 'ftp://optics.gi.alaska.edu'


def download(startend: Tuple[datetime, datetime],
             site: str, odir: Path,
             host: str = None,
             wavelen: str = None,
             overwrite: bool = False) -> List[Path]:
    """
    startend: tuple of datetime
    year,month,day: integer
    hour, minute:  start,stop integer len == 2
    """
    if not host:
        host = HOST

    assert len(startend) == 2

    start = parse(startend[0]) if isinstance(startend[0], str) else startend[0]  # type: ignore
    end = parse(startend[1]) if isinstance(startend[1], str) else startend[1]  # type: ignore

    if end < start:
        raise ValueError('start time must be before end time!')

    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2] + site
    odir = Path(odir).expanduser().resolve()
    odir.mkdir(exist_ok=True)
# %% get available files for this day
    rparent = f'{fpath}/DASC/RAW/{start.year:4d}'
    rday = f'{start.year:4d}{start.month:02d}{start.day:02d}'
# %% wavelength
    if wavelen is None:
        pass
    elif isinstance(wavelen, int):
        wavelen = f'{wavelen:04d}'
    elif isinstance(wavelen, str):
        if len(wavelen) != 4:
            raise ValueError('expecting 4-character wavelength spec e.g. 0428')
    elif not isinstance(wavelen, (tuple, list)):
        raise TypeError('expecting 4-character wavelength spec e.g. 0428')

    flist = []

    with ftplib.FTP(ftop, 'anonymous', 'guest', timeout=15) as F:
        F.cwd(rparent)
        dlist = F.nlst()
        if rday not in dlist:
            raise FileNotFoundError(f'{rday} does not exist under {host}/{rparent}')

        print('downloading to', odir)
        F.cwd(rday)
        dlist = F.nlst()

        print(f'remote filesize approx. {F.size(dlist[0])/1000} kB.')  # type: ignore

        for rfn in dlist:
            # %% qualifiers
            if wavelen and rfn[9:13] not in wavelen:
                continue

            tfile = datetime.strptime(rfn[14:-9], '%Y%m%d_%H%M%S')
            if tfile < start or tfile > end:
                continue
# %% download file
            ofn = odir / rfn
            flist.append(ofn)

            if not overwrite:
                if ofn.is_file():  # do filesizes match, if so, skip download
                    rsize = F.size(rfn)
                    if ofn.stat().st_size == rsize:
                        print('SKIPPING existing', ofn)
                        continue

            print(ofn)
            with ofn.open('wb') as h:
                F.retrbinary(f'RETR {rfn}', h.write)
                sleep(1)  # anti-leech

    return flist
