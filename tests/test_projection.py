#!/usr/bin/env python
from pathlib import Path
import dascutils as du
import pytest
from pytest import approx

R = Path(__file__).parent
azelstem = R.parent/'cal/PKR_DASC_20110112'


def test_projection():
    dp = pytest.importorskip('dascutils.projection')

    data = du.load(R, azelstem)

    data = dp.project_altitude(data, 100.)

    assert data.mapping_alt_km == approx(100.)
    assert data.mapping_lat[266, 247] == approx(65.12351)
    assert data.mapping_lon[266, 247] == approx(-147.48196)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
