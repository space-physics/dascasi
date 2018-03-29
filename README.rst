.. image:: https://zenodo.org/badge/51016067.svg
   :target: https://zenodo.org/badge/latestdoi/51016067

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

.. image:: https://travis-ci.org/scivision/dascutils.svg?branch=master
    :target: https://travis-ci.org/scivision/dascutils

.. image:: https://coveralls.io/repos/github/scivision/dascutils/badge.svg?branch=master
    :target: https://coveralls.io/github/scivision/dascutils?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/xrtb6fc3d4ojp507?svg=true
    :target: https://ci.appveyor.com/project/scivision/dascutils

.. image:: https://api.codeclimate.com/v1/badges/36b08deedc7d2bf750c8/maintainability
   :target: https://codeclimate.com/github/scivision/dascutils/maintainability
   :alt: Maintainability

============
DASC utils
============

Utilities for plotting, saving, analyzing the Poker Flat Research Range Digital All Sky Camera.
(Other locations, too).

This program handles the corrupted FITS files due to the RAID array failure on 2013 data.

The raw data FITS are one image per file.

.. contents::

Install
=======
::

	pip install -e .

Download raw DASC files by time
===========================
Example download October 7, 2015 from 8:23 to 8:54 UTC::

    ./DownloadDASC.py 2015-10-07T08:23Z 2015-10-07T08:54Z

-o  download directory
-c  clobber existing files
-s  three-letter site acronym PKR for poker flat etc.

Make movies from DASC raw data files
====================================
Plots all wavelengths in subplots::

    ./PlotDASC.py datadir

Spatial registration (plate scale)
==================================
the ``cal/`` directory contains ``AZ`` and ``EL`` files corresponding to each pixel.
Be sure you know if you're using magnetic north or geographic north, or you'll see a rotation by the declination.

Note the date in the filename--perhaps the camera was moved since before or long after that date?
