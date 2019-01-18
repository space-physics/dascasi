#!/usr/bin/env python
from pathlib import Path
import xarray
import pytest
from datetime import datetime
#
import dascutils as du

R = Path(__file__).parent

azelstem = R.parent/'cal/PKR_DASC_20110112'


def test_nonexistent_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        du.load(tmp_path)


@pytest.mark.parametrize('wavelength, L, t',
                         [('0428', 1, datetime(2015, 10, 7, 8, 23, 55, 961000)),
                          ('0558', 2, datetime(2015, 10, 7, 8, 23, 51, 743000)),
                          ('0630', 1, datetime(2015, 10, 7, 8, 23, 59, 586000))])
def test_basic_load(wavelength, L, t):
    dset = du.load(R)
    assert isinstance(dset, xarray.Dataset)

    data = dset[wavelength].dropna(dim='time')
    assert data.shape == (L, 512, 512)

    assert 'az' not in dset.data_vars
    assert dset.lat == 65.126
    assert dset.lon == -147.479

    assert data.time.values[0].astype('datetime64[us]').astype(datetime) == t


def test_timerange_and_wavelength():
    data = du.load(R, azelstem,
                   treq=('2012-01-03T08:32:02', '2016-01-04'),
                   wavelenreq='0558')
    assert data['0558'].shape == (2, 512, 512)


@pytest.mark.parametrize('wavelength, L', [('0558', 1)])
def test_singletime(wavelength, L):
    data = du.load(R, azelstem, '2012-01-03T08:32:02')
    assert data[wavelength].shape == (L, 512, 512)
    assert data.az.shape == (512, 512)
    assert data.el.shape == (512, 512)
    assert data.time.values.astype('datetime64[us]').astype(datetime) == datetime(2015, 10, 7, 8, 23, 51, 743000)


@pytest.mark.parametrize('wavelength, L', [('0428', 1), ('0558', 2), ('0630', 1)])
def test_full_load(wavelength, L):
    # %% wavelength request
    data = du.load(R, azelstem, wavelenreq=wavelength)
    assert data[wavelength].shape == (L, 512, 512)


def test_write_hdf5(tmp_path):
    pytest.importorskip('netCDF4')

    ofn = tmp_path / 'test.nc'

    ref = du.load(R, azelstem, ofn=ofn)

    dat = xarray.open_dataset(ofn)

    assert dat.equals(ref)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
