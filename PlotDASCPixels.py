#!/usr/bin/env python
"""
Plots time series of pixel brightness

example plot time series of pixel brightness projected to 100 km at lat, lon 65.3 -148.2
python PlotDASCPixels.py tests/ 100. 65.3 -148.2
"""

from pathlib import Path
from matplotlib.pyplot import show
import dascutils as du

R = Path(__file__).parent
data_dir = R / "tests"
cal_stem = R / "cal" / "PKR_DASC_20110112_"

mapping_altitude_km = {"0428": 110.0, "0558": 150.0, "0630": 200.0, "unknown": 150.0}

imgs = du.load(data_dir, cal_stem, wavelength_altitude_km=mapping_altitude_km)

show()
