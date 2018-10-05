import pytest
import dascutils as du
import tempfile
import socket
import subprocess
from pathlib import Path

R = Path(__file__).parent


def test_mod():

    try:
        with tempfile.TemporaryDirectory() as d:
            flist = du.download(('2015-10-07T08:23:54', '2015-10-07T08:23:56'),
                                site='PKR', odir=d)
            for fn in flist:
                assert fn.is_file()
                assert fn.parent.samefile(d)

    except socket.gaierror as e:
        pytest.skip(f"Bad internet connection?   {e}")


def test_script():
    try:
        subprocess.check_call(['DownloadDASC',
                               '2015-10-07T08:23:54',
                               '2015-10-07T08:23:56', str(R)])
    except socket.gaierror as e:
        pytest.skip(f"Bad internet connection?   {e}")


if __name__ == '__main__':
    pytest.main(['-x', __file__])
