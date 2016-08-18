.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/
.. image:: https://travis-ci.org/scienceopen/dascutils.svg?branch=master
    :target: https://travis-ci.org/scienceopen/dascutils

.. image:: https://coveralls.io/repos/github/scienceopen/dascutils/badge.svg?branch=master 
    :target: https://coveralls.io/github/scienceopen/dascutils?branch=master    

============
DASC utils
============

Utilities for plotting, saving, analyzing the Poker Flat Research Range Digital All Sky Camera.
(Other locations, too).

This program handles the corrupted FITS files due to the RAID array failure on 2013 data.

The raw data FITS are one image per file

.. contents::

Install
=======
::

	python setup.py develop
	
Download raw DASC files by time
===========================
Example download October 7, 2015 from 8:23 to 8:54 UTC::

    ./DownloadDASC.py 2015 10 7 8:23 8:54 
    
-o  download directory
-c  clobber existing files
-s  three-letter site acronym PKR for poker flat etc.

Make movies from DASC raw data files
====================================
Plots all wavelengths in subplots::

    ./PlotDASC.py datadir
