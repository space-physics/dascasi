#!/usr/bin/env python
from pathlib import Path
import dascutils as du
import dascutils.readDASCfits
import unittest
import numpy as np

R = Path(__file__).parent
# .resolve() was added for Travis-CI debugging, not necessary
fn = (R/'PKR_DASC_0428_20151007_082305.930.FITS').resolve()
azfn = (R.parent/'cal/PKR_DASC_20110112_AZ_10deg.fits').resolve()
elfn = (R.parent/'cal/PKR_DASC_20110112_EL_10deg.fits').resolve()

assert fn.is_file(),f'could not find test data file {fn}'
assert azfn.is_file(),f'could not find azimuth cal file {azfn}'
assert elfn.is_file(),f'could not find elevation cal file {azfn}'


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


    def test_getdasc(self):
        du.getdasc(('2015-10-07T08:23:04','2015-10-07T08:23:06'),
                'ftp://optics.gi.alaska.edu', 'PKR', R)



if __name__ == '__main__':
    unittest.main()