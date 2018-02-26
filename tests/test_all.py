#!/usr/bin/env python
from pathlib import Path
import dascutils as du
import dascutils.readDASCfits
import unittest
import numpy as np

R = Path(__file__).parent
fn = R/'PKR_DASC_0428_20151007_082305.930.FITS'
azfn = R.parent/'cal/PKR_DASC_20110112_AZ_10deg.fits'
elfn = azfn = R.parent/'cal/PKR_DASC_20110112_EL_10deg.fits'


class BasicTest(unittest.TestCase):

    def test_timestamp(self):
        np.testing.assert_allclose(du.totimestamp('2012-01-03T08:32:02Z'), 1325579522.0)
        np.testing.assert_allclose(du.totimestamp(1325579522), 1325579522.0)
        np.testing.assert_allclose(du.totimestamp(['2012-01-03T08:32:02Z','2012-01-03T08:32:12Z']),
                                   (1325579522.0, 1325579532.0))

    def test_readdasc(self):

        img,times,az,el,sensorloc,wlused = du.readDASCfits.readallDasc(fn,azfn,elfn,None)
        assert isinstance(img[0], np.ndarray)
        assert img[0].shape == (1,512,512)
        assert sensorloc == {'lat': 65.126, 'lon': -147.479, 'alt_m': 200.0}

if __name__ == '__main__':
    unittest.main()
