import pytest
import dascutils as du
import socket
import ftplib
from pathlib import Path

R = Path(__file__).parent


def test_nonexist_site(tmp_path):
    with pytest.raises(ValueError):
        du.download(
            ("2015-10-07T08:23:50", "2015-10-07T08:23:56"), site="pk", odir=tmp_path, wavelen="428"
        )


def test_nonexist_remote(tmp_path):
    with pytest.raises(ValueError):
        du.download(
            ("2015-10-07T08:23:50", "2015-10-07T08:23:56"), site="PKR", odir=tmp_path, wavelen="428"
        )


@pytest.mark.parametrize(
    "wavelength", (["0428"], ("0428", "0558")), ids=("one_wavelength", "two_wavelengths")
)
def test_mod(tmp_path, wavelength):

    dpath = tmp_path

    try:

        flist = du.download(
            ("2015-10-07T08:23:50", "2015-10-07T08:23:56"),
            site="PKR",
            odir=dpath,
            wavelen=wavelength,
        )
        for fn in flist:
            assert fn.is_file()
            assert fn.parent.samefile(dpath)

        assert len(flist) == len(wavelength)

    except (socket.gaierror, socket.timeout, ftplib.error_temp) as e:
        pytest.skip(f"Bad internet connection?   {e}")
