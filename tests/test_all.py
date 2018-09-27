#!/usr/bin/env python
from pathlib import Path
import numpy as np
import xarray
import tempfile
import pytest
from pytest import approx
from datetime import datetime
#
import dascutils as du
import dascutils.io as dio

R = Path(__file__).parent

azelstem = R.parent/'cal/PKR_DASC_20110112'


def test_timestamp():
    assert du.totimestamp('2012-01-03T08:32:02Z') == approx(1325579522.0)
    assert du.totimestamp(1325579522) == approx(1325579522.0)
    assert du.totimestamp(['2012-01-03T08:32:02Z', '2012-01-03T08:32:12Z']) == approx((1325579522.0, 1325579532.0))


def test_basic_load():
    # %% expected fail
    with pytest.raises(FileNotFoundError):
        with tempfile.TemporaryDirectory() as d:
            data = dio.load(d)
# %% most basic
    data = dio.load(R)
    assert isinstance(data, xarray.Dataset)

    d428 = data[428].dropna(dim='time')
    assert d428.shape == (1, 512, 512)

    d558 = data[558].dropna(dim='time')
    assert d558.shape == (2, 512, 512)

    d630 = data[630].dropna(dim='time')
    assert d630.shape == (1, 512, 512)

    assert 'az' not in data.data_vars
    assert data.lat == 65.126
    assert data.lon == -147.479

    assert d630.time.values == np.datetime64('2015-10-07T08:23:59.586000')
    assert data.cadence[558].values == np.timedelta64(12500, 'ms')


def test_full_load():
    # %% single time request
    data = dio.load(R, azelstem, '2012-01-03T08:32:02')
    assert data[558].shape == (1, 512, 512)
    assert data.az.shape == (512, 512)
    assert data.el.shape == (512, 512)
    assert data.time.item() == datetime(2015, 10, 7, 8, 23, 51, 743000)
# %% multi-time request and wavelength
    data = dio.load(R, azelstem, ('2012-01-03T08:32:02', '2016-01-04'), 558)
    assert data[558].shape == (2, 512, 512)
# %% wavelength request
    data = dio.load(R, azelstem, None, 428)
    assert data[428].shape == (1, 512, 512)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
