#!/usr/bin/env python
"""
convert DASC FITS stack to HDF5
"""
import dascutils as du
from pathlib import Path
import argparse


p = argparse.ArgumentParser()
p.add_argument("inpath", help="Directory where DASC FITS files are")
p.add_argument("outfile", help="filename to write .h5 HDF5 file to")
p.add_argument("-t", "--tlim", help="start stop time bounds", nargs=2)
P = p.parse_args()


outfile = Path(P.outfile).expanduser()
if outfile.is_dir():
    raise IsADirectoryError(outfile)
outfile.parent.mkdir(exist_ok=True, parents=True)

imgs = du.load(P.inpath, treq=P.tlim)
du.save_hdf5(imgs, outfile)
