#!/usr/bin/env python
from pathlib import Path
import numpy as np
import xarray
import tempfile
from datetime import datetime
from numpy.testing import assert_allclose
#
import dascutils as du
import dascutils.io as dio

R = Path(__file__).parent

fn = R/'PKR_DASC_0428_20151007_082305.930.FITS'
azfn = R.parent/'cal/PKR_DASC_20110112_AZ_10deg.fits'
elfn = R.parent/'cal/PKR_DASC_20110112_EL_10deg.fits'

assert fn.is_file(),f'could not find test data file {fn}'
assert azfn.is_file(),f'could not find azimuth cal file {azfn}'
assert elfn.is_file(),f'could not find elevation cal file {azfn}'


def test_timestamp():
    assert_allclose(du.totimestamp('2012-01-03T08:32:02Z'), 1325579522.0)
    assert_allclose(du.totimestamp(1325579522), 1325579522.0)
    assert_allclose(du.totimestamp(['2012-01-03T08:32:02Z','2012-01-03T08:32:12Z']),
                                   (1325579522.0, 1325579532.0))

def test_lost():
    data = dio.load(fn,azfn,elfn)
    assert isinstance(data, xarray.Dataset)
    assert 428 in data.data_vars
    assert data[428].shape == (1,512,512)

    assert data.lat == 65.126
    assert data.lon == -147.479
    assert data.time.values == datetime(2015, 10, 7, 8, 23, 5, 930000)
# %% no time request
    data = dio.load(fn)
    assert data[428].shape == (1,512,512)
    assert 'az' not in data.data_vars
# %% single time request
    data = dio.load(fn,azfn,elfn,'2012-01-03T08:32:02')
    assert data[428].shape == (1,512,512)
    assert data.az.shape == (512,512)
    assert data.el.shape == (512,512)
# %% multi-time request
    data = dio.load(fn,azfn,elfn,('2012-01-03T08:32:02','2016-01-04'))
    assert data[428].shape == (1,512,512)
# %% wavelength request
    data = dio.load(fn,azfn,elfn,None,428)
    assert data[428].shape == (1,512,512)


def atest_download():
    with tempfile.TemporaryDirectory() as d:
        du.download(('2015-10-07T08:23:04','2015-10-07T08:23:06'),
                    'ftp://optics.gi.alaska.edu', 'PKR', d)



if __name__ == '__main__':
    np.testing.run_module_suite()
