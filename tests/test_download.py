import pytest
import dascutils as du
from pathlib import Path

R = Path(__file__).parent

def test_download():

    try:
        du.download(('2015-10-07T08:23:54', '2015-10-07T08:23:56'),
                    'ftp://optics.gi.alaska.edu', 'PKR', R)
    except Exception as e:
        pytest.skip(f"Bad internet connection?   {e}")
        

if __name__ == '__main__':
    pytest.main(['-x', __file__])
