#!/usr/bin/env python
from setuptools import setup

req = ['nose','python-dateutil','pytz','numpy','astropy','scipy','matplotlib',
        'sciencedates','themisasi',]

setup(name='dascutils',
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/dascutils',
      description='Utilities for UAF Digital All-Sky Camera: reading and plotting',
      version = '0.5',
      classifiers=[
      'Programming Language :: Python :: 3.6',
      ],
      install_requires=req,
     dependency_links = [
        'https://github.com/scivision/themisasi/tarball/master#egg=themisasi-999.0.0'
	],
#      setup_requires=['numpy'], #for spacepy/themisasi doesn't work
      packages=['dascutils'],
	  )


