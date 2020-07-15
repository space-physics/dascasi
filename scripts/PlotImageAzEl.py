#!/usr/bin/env python
"""
Plots image with azimuth/elevation overlay

for simplicity, hard-coded to test data, but can be used for any data.
"""

from pathlib import Path
from matplotlib.pyplot import show
import dascutils as du
import dascutils.plots as dup

R = Path(du.__file__).parent
data_dir = R / "tests/data/"
data_file = data_dir / "PKR_DASC_0558_20151007_082351.743.FITS"
cal_stem = data_dir / "cal/PKR_DASC_0558_20150213_"


imgs = du.load(data_file, cal_stem)

dup.image_azel(imgs, cal_stem.name, data_file.name)

show()
