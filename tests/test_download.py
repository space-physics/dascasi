import pytest
import dascutils as du
import tempfile
import socket
import ftplib
import subprocess
from pathlib import Path

R = Path(__file__).parent


def test_mod():

    with pytest.raises(ValueError):
        du.download(('2015-10-07T08:23:50', '2015-10-07T08:23:56'),
                    site='PKR', odir=R, wavelen='428')

    try:
        with tempfile.TemporaryDirectory() as d:

            flist = du.download(('2015-10-07T08:23:50', '2015-10-07T08:23:56'),
                                site='PKR', odir=d, wavelen='0428')
            for fn in flist:
                assert fn.is_file()
                assert fn.parent.samefile(d)

            assert len(flist) == 1

            flist = du.download(('2015-10-07T08:23:50', '2015-10-07T08:23:56'),
                                site='PKR', odir=d, wavelen=['0428', '0558'])
            for fn in flist:
                assert fn.is_file()
                assert fn.parent.samefile(d)

            assert len(flist) == 2

    except (socket.gaierror, ftplib.error_temp) as e:
        pytest.skip(f"Bad internet connection?   {e}")


def test_script():
    try:
        subprocess.check_call(['DownloadDASC',
                               '2015-10-07T08:23:54',
                               '2015-10-07T08:23:56', str(R)])
    except subprocess.CalledProcessError as e:
        pytest.skip(f"Bad internet connection?   {e}")


if __name__ == '__main__':
    pytest.main(['-x', __file__])
