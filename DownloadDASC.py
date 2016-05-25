#!/usr/bin/env python
import ftplib
from urllib.parse import urlparse

host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"

def getdasc(path=''):
    parsed = urlparse(host)
    ftop = parsed[1]
    fpath = parsed[2]

    with ftplib.FTP(ftop,'anonymous','guest',timeout=5) as F:
        F.cwd(fpath)
        flist = F.nlst()
        print(flist)
    # TODO filter times, FTP site was down at this time

if __name__ == '__main__':
    getdasc()