#!/usr/bin/env python
"""
Plots image projected to per-wavelength altitude

for simplicity, hard-coded to test data, but can be used for any data.
"""

from pathlib import Path
from matplotlib.pyplot import show
import dascutils as du
import dascutils.plots as dup

R = Path(du.__file__).parent
data_dir = R / "tests/data"
cal_stem = R / "cal/PKR_DASC_20110112_"


mapping_altitude_km = {"0428": 110.0, "0558": 150.0, "0630": 200.0, "0000": 150.0}

imgs = du.load(data_dir, cal_stem, wavelength_altitude_km=mapping_altitude_km)

for wavelength in mapping_altitude_km.keys():
    if wavelength not in imgs:
        continue
    dup.plot_projected_image(imgs[wavelength])

show()
