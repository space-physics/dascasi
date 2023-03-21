#!/usr/bin/env python
"""
Plots image with azimuth/elevation overlay

for simplicity, hard-coded to test data, but can be used for any data.
"""

from pathlib import Path
from matplotlib.pyplot import show

import dascasi as du
import dascasi.plots as dup

R = Path(du.__file__).parent
data_dir = R / "tests/data/"

# data_file = data_dir / "PKR_DASC_0558_20151007_082351.743.FITS"

# %% full moon
# dascasi_download.exe PKR 2015-10-27T12:23:19 2015-10-27T12:23:21 ~/data/dasc
data_file = Path("~/data/dasc/PKR_DASC_0428_20151027_122319.645.FITS")

cal_stem = data_dir / "cal/PKR_DASC_0558_20150213_"


imgs = du.load(data_file, cal_stem)

# rearrange data to fit physical world
imgs["az"] = 360 - imgs["az"]  # empirical
assert imgs["0428"].shape == (1, 512, 512)
for k in ["az", "el"] + imgs["wavelengths"].tolist():
    if imgs[k].ndim == 2:
        imgs[k] = imgs[k].transpose(1, 0)  # empirical
    else:
        imgs[k] = imgs[k].transpose("time", "x", "y")  # empirical

dup.image_azel(imgs, cal_stem.name, data_file.name)

show()
