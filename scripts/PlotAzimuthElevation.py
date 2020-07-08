#!/usr/bin/env python
"""
To get azimuth, elevation per pixel, include the optional second argument to du.load().
The azimuth, elevation data is not validated or checked.
If the calibration data is from a different camera or camera position, the azimuth, elevation data will be incorrect.

The quirk about specifying the filename cal_stem is that azimuth,elevation data are given as two files.
The stem is the common part of the fullpath to these files.
"""

from pathlib import Path
import dascutils as du
import dascutils.plots as dup
from matplotlib.pyplot import show

R = Path(__file__).parent

data_dir = R / "../src/dascutils/tests/data"
cal_stem = R / "../src/dascutils/cal/PKR_DASC_20110112"

imgs = du.load(data_dir, cal_stem)
# If you omit "cal_stem", `imgs` will not have "az,el" keys

dup.contour_azel(imgs, cal_stem)
show()
