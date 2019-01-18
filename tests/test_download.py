import pytest
import dascutils as du
import socket
import ftplib
import subprocess
from pathlib import Path

R = Path(__file__).parent


def test_nonexistent_remote():
    with pytest.raises(ValueError):
        du.download(('2015-10-07T08:23:50', '2015-10-07T08:23:56'),
                    site='PKR', odir=R, wavelen='428')


@pytest.mark.parametrize('wavelength',
                         (['0428'], ('0428', '0558')),
                         ids=('one_wavelength', 'two_wavelengths'))
def test_mod(tmp_path, wavelength):

    dpath = tmp_path

    try:

        flist = du.download(('2015-10-07T08:23:50', '2015-10-07T08:23:56'),
                            site='PKR', odir=dpath, wavelen=wavelength)
        for fn in flist:
            assert fn.is_file()
            assert fn.parent.samefile(dpath)

        assert len(flist) == len(wavelength)

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
