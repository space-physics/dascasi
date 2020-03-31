#!/usr/bin/env python
from pathlib import Path
import xarray
import pytest
import numpy as np
from pytest import approx
from datetime import datetime

import dascutils as du

R = Path(__file__).parent

azelstem = R.parent / "cal/PKR_DASC_20110112"


def test_nonexistent_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        du.load(tmp_path)


@pytest.mark.parametrize(
    "wavelength, L, t",
    [
        ("0428", 1, datetime(2015, 10, 7, 8, 23, 55, 961000)),
        ("0558", 2, datetime(2015, 10, 7, 8, 23, 51, 743000)),
        ("0630", 1, datetime(2015, 10, 7, 8, 23, 59, 586000)),
    ],
)
def test_basic_load(wavelength, L, t):
    imgs = du.load(R / "data")
    assert isinstance(imgs, dict)
    assert isinstance(imgs[wavelength], xarray.DataArray)

    data = imgs[wavelength]
    assert data.shape == (L, 512, 512)
    assert imgs["lat0"] == approx(65.126)
    assert imgs["lon0"] == approx(-147.479)

    assert data.time.values[0].astype("datetime64[us]").astype(datetime) == t


def test_timerange_and_wavelength():
    data = du.load(R / "data", azelstem, treq=("2012-01-03T08:32:02", "2016-01-04"), wavelenreq="0558")
    assert data["0558"].shape == (2, 512, 512)
    assert "0428" not in data
    assert "0630" not in data


@pytest.mark.parametrize("wavelength, L", [("0558", 1)])
def test_singletime(wavelength, L):
    data = du.load(R / "data", azelstem, treq="2012-01-03T08:32:02")
    assert data[wavelength].shape == (L, 512, 512)
    assert data["az"].shape == (512, 512)
    assert data["el"].shape == (512, 512)
    assert data[wavelength].time.values.astype("datetime64[us]").astype(datetime) == datetime(2015, 10, 7, 8, 23, 51, 743000)


@pytest.mark.parametrize("wavelength, L", [("0428", 1), ("0558", 2), ("0630", 1)])
def test_full_load(wavelength, L):
    # %% wavelength request
    data = du.load(R / "data", azelstem, wavelenreq=wavelength)
    assert data[wavelength].shape == (L, 512, 512)


def test_read_write_hdf5(tmp_path):
    outfn = tmp_path / "test.h5"

    ref = du.load(R / "data", azelstem)
    du.save_hdf5(ref, outfn)

    dat = du.load(outfn)

    for k, v in dat.items():
        if isinstance(v, np.ndarray):
            assert (v == ref[k]).all()
        elif isinstance(v, xarray.DataArray):
            assert v.equals(ref[k])
        else:
            assert v == ref[k]


if __name__ == "__main__":
    pytest.main([__file__])
