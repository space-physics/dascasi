#!/usr/bin/env python
from pathlib import Path
import xarray
import tempfile
import pytest
from pytest import approx
from datetime import datetime
#
import dascutils as du

R = Path(__file__).parent

azelstem = R.parent/'cal/PKR_DASC_20110112'


def test_basic_load():
    # %% expected fail
    with pytest.raises(FileNotFoundError):
        with tempfile.TemporaryDirectory() as d:
            data = du.load(d)
# %% most basic
    data = du.load(R)
    assert isinstance(data, xarray.Dataset)

    d428 = data['0428'].dropna(dim='time')
    assert d428.shape == (1, 512, 512)

    d558 = data['0558'].dropna(dim='time')
    assert d558.shape == (2, 512, 512)

    d630 = data['0630'].dropna(dim='time')
    assert d630.shape == (1, 512, 512)

    assert 'az' not in data.data_vars
    assert data.lat == 65.126
    assert data.lon == -147.479

    assert d630.time.values.astype('datetime64[us]').astype(datetime) == datetime(2015, 10, 7, 8, 23, 59, 586000)


def test_full_load():
    # %% single time request
    data = du.load(R, azelstem, '2012-01-03T08:32:02')
    assert data['0558'].shape == (1, 512, 512)
    assert data.az.shape == (512, 512)
    assert data.el.shape == (512, 512)
    assert data.time.values.astype('datetime64[us]').astype(datetime) == datetime(2015, 10, 7, 8, 23, 51, 743000)
# %% multi-time request and wavelength
    data = du.load(R, azelstem,
                   treq=('2012-01-03T08:32:02', '2016-01-04'),
                   wavelenreq='0558')
    assert data['0558'].shape == (2, 512, 512)
# %% wavelength request
    data = du.load(R, azelstem, wavelenreq='0428')
    assert data['0428'].shape == (1, 512, 512)


def test_projection():
    pytest.importorskip('pymap3d')
    pytest.importorskip('scipy')

    data = du.load(R, azelstem)

    data = du.project_altitude(data, 100.)

    assert data.mapping_alt_km == approx(100.)
    assert data.mapping_lat[266, 247] == approx(65.12351)
    assert data.mapping_lon[266, 247] == approx(-147.48196)


def test_write_hdf5():
    pytest.importorskip('netcdf4')

    with tempfile.TemporaryDirectory() as D:

        ofn = Path(D) / 'test.nc'

        ref = du.load(R, azelstem, ofn=ofn)

        dat = xarray.open_dataset(ofn, autoclose=True)

        assert dat.equals(ref)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
