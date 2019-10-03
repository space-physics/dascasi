[![image](https://zenodo.org/badge/51016067.svg)](https://zenodo.org/badge/latestdoi/51016067)

[![Actions Status](https://github.com/space-physics/dascutils/workflows/ci_python/badge.svg)](https://github.com/space-physics/dascutils/actions)

[![PyPi versions](https://img.shields.io/pypi/pyversions/dascutils.svg)](https://pypi.python.org/pypi/dascutils)
[![PyPi Download stats](http://pepy.tech/badge/dascutils)](http://pepy.tech/project/dascutils)


# DASC all-sky camera utilities

Utilities for plotting, saving, analyzing the Poker Flat Research Range Digital All Sky Camera. (Other locations, too).

This program handles the corrupted FITS files due to the RAID array failure on 2013 data.

The raw data FITS are one image per file.


## Install

Most people will find it useful to have the example scripts and the tests built into the Git repo.

```sh
git clone https://github.com/space-physics/dascutils/actions

pip install -e dascutils
```

## Usage
Many analysts may use the API directly, like:
```python
import dascutils as du

data = du.load('tests/PKR_DASC_0558_20151007_082351.743.FITS')
```
This returns an [xarray.Dataset](http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html), which is like a "smart" Numpy array.
The images are index by wavelength if it was specified in the data file, or 'unknown' otherwise.
The images are in a 3-D stack: (time, x, y).
`data.time` is the time of each image.
also several metadata parameters are included like the location of the camera.

Save the data using lossless compression to NetCDF4 / HDF5 by

```python
du.save_hdf5(data, "foo.nc")
```

Now we give several examples.

### Download raw DASC files by time

Download Poker Flat Research Range "PKR" October 7, 2015 from 8:23 to 8:54 UTC to `~/data/`:

```sh
python dascutils/DownloadDASC.py PKR 2015-10-07T08:23 2015-10-07T08:54 ~/data
```

* `-w` four-letter wavelength in nanometers e.g. 0630

As usual, we assume UTC and do NOT specify the timezone.

### convert FITS stack to HDF5 / NetCDF4

It is very tedious to download large amounts of DASC data in single FITS files.
We have tried to make this faster by multi-threading the download, but then the FTP server anti-leeching
leaves us with broken downloads.
As an alternative in general, it's more convenient to have a single HDF5 file for a day rather than 10,000 FITS files.
Convert a bunch of FITS files to HDF5 like:

```sh
python dascutils/ConvertDASC_FITS_to_HDF5.py ~/data/2015-10-07 ~/data/2015-10-07.nc
```

* `-t` start stop times to convert

### Make movies from DASC raw data files

Play movie of all wavelengths in subplots for files in a directory, for example:

```sh
python dascutils/PlayMovie.py dascutils/tests/ -a cal/PKR_DASC_20110112
```

additional options include:

* `-t` specifiy time limits e.g.  `-t 2014-01-02T02:30 2014-01-02T02:35`
* `-w` choose only certain wavelength(s)

### Spatial registration (plate scale)

The `cal/` directory contains `AZ` and `EL` files corresponding to each pixel.

```python
import dascutils as du

data = du.load('tests/PKR_DASC_0558_20151007_082351.743.FITS', azelfn='cal/PKR_DASC_20110112')
```

now `data` includes data variables `az` and `el`, same shape as the image(s), along with camera position in `lat` `lon` `alt_m`.

* Be sure you know if you're using magnetic north or geographic north, or you'll see a rotation by the declination.
* Note the date in the filename--perhaps the camera was moved since before or long after that date?

### Map Projection
A common task in auroral and airglow analysis is to project the image to an imaginary alttiude, that is, as if all the brightness were coming from that altitude.
Typically that altitude is on the order of 100 km.
The `dascutils.project_altitude()` function adds coordinates `mapping_lat` `mapping_lon` to the xarray.Dataset by:

```python
import dascutils as du
import dascutils.projection as dp

data = du.load('myfile.FITS', azelfn='cal/PKR_DASC_20110112')

data = dp.project_altitude(data, 100.)  # for 100 km
```

The `dascutils.projection` is a separate import because it calls extra Python modules that aren't needed for basic data loading.

