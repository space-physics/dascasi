from .io import load

from pathlib import Path
from argparse import ArgumentParser


def _plotdasc(img, outdir: Path, cadence: float):
    from .plots import histogram_dasc, moviedasc

    histogram_dasc(img, outdir)

    moviedasc(img, outdir, cadence)


def dascasi_movie():
    from matplotlib.pyplot import show

    p = ArgumentParser(
        description="for DASC all sky camera, read and play az/el mapping and images"
    )
    p.add_argument("indir", help="directory of .fits or specific .fits file")
    p.add_argument("-o", "--outdir", help="directory to write plots to")
    p.add_argument("-t", "--tlim", help="only plot data in this range", nargs=2)
    p.add_argument("-a", "--azelfn", help="stem for DASC .fits azimuth calibration")
    p.add_argument(
        "-w",
        "--wavelength",
        help="select wavelength(s) to plot simultaneously [428 558 630]",
        nargs="+",
    )
    p.add_argument(
        "-c",
        "--cadence",
        help="set playback cadence to request times [sec]",
        type=float,
        default=5.0,
    )
    p = p.parse_args()

    imgs = load(p.indir, p.azelfn, treq=p.tlim, wavelenreq=p.wavelength)

    _plotdasc(imgs, p.outdir, p.cadence)

    show()


if __name__ == "__main__":
    dascasi_movie()
