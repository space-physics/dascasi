import argparse

from .web import download


def dascasi_download():
    """
    Download DASC ASI data for various sites.

    site VEE for one minute:

        python DownloadDASC.py VEE 2017-02-13T06:59 2017-02-13T07:00 ~/data
    """
    p = argparse.ArgumentParser(description="download DASC all-sky camera data")
    p.add_argument("site", choices=["EAA", "FYU", "KAK", "PKR", "TOO", "VEE"])
    p.add_argument(
        "startend", help="start/end times UTC e.g. 2012-11-03T06:23 2012-11-03T07", nargs=2
    )
    p.add_argument("odir", help="directory to write downloaded FITS to")
    p.add_argument("-w", "--wavelen", help="request specific wavelength(s)", nargs="+")
    p.add_argument("-host", default="ftp://optics.gi.alaska.edu")
    p = p.parse_args()

    # host = "ftp://mirrors.arsc.edu/AMISR/PKR/DASC/RAW/"
    download(p.startend, p.site, p.odir, p.host, p.wavelen)


if __name__ == "__main__":
    dascasi_download()
