#!/usr/bin/env python
from pathlib import Path
import dascutils as du
import pytest
from pytest import approx

R = Path(__file__).parent
azelstem = R / "data/cal/PKR_DASC_0558_20150213_"
mapping_altitude_km = {"0428": 110.0, "0558": 150.0, "0630": 200.0}


@pytest.mark.parametrize(
    "wavelen,refalt,reflat,reflon",
    [
        ("0428", 110.0, 65.34932617, -147.92703665),
        ("0558", 150.0, 65.42698381, -148.0819881),
        ("0630", 200.0, 65.52119085, -148.29018938),
    ],
)
def test_projection(wavelen, refalt, reflat, reflon):

    data = du.load(R / "data", azelstem, wavelenreq=wavelen, wavelength_altitude_km=mapping_altitude_km)

    assert data[wavelen].mapping_alt_km == approx(refalt)

    assert data[wavelen].lat[266, 247] == approx(reflat)
    assert data[wavelen].lon[266, 247] == approx(reflon)


if __name__ == "__main__":
    pytest.main([__file__])
