[![image](https://zenodo.org/badge/51016067.svg)](https://zenodo.org/badge/latestdoi/51016067)
[![image](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![image](https://travis-ci.org/scivision/dascutils.svg?branch=master)](https://travis-ci.org/scivision/dascutils)
[![image](https://coveralls.io/repos/github/scivision/dascutils/badge.svg?branch=master)](https://coveralls.io/github/scivision/dascutils?branch=master)
[![image](https://ci.appveyor.com/api/projects/status/xrtb6fc3d4ojp507?svg=true)](https://ci.appveyor.com/project/scivision/dascutils)
[![Maintainability](https://api.codeclimate.com/v1/badges/36b08deedc7d2bf750c8/maintainability)](https://codeclimate.com/github/scivision/dascutils/maintainability)

# DASC all-sky camera utilitiess

Utilities for plotting, saving, analyzing the Poker Flat Research Range Digital All Sky Camera. (Other locations, too).

This program handles the corrupted FITS files due to the RAID array failure on 2013 data.

The raw data FITS are one image per file.


## Install

```sh
pip install -e .
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

### Download raw DASC files by time

Example download October 7, 2015 from 8:23 to 8:54 UTC to `~/data/`:

```sh
DownloadDASC 2015-10-07T08:23 2015-10-07T08:54 ~/data
```

* `-c` overwrite existing files
* `-s` three-letter site acronym e.g. `PKR` for poker flat etc.

### Make movies from DASC raw data files

Plots all wavelengths in subplots, for example:

```sh
PlotDASC tests/ -a cal/PKR_DASC_20110112
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

