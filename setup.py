#!/usr/bin/env python
req = ['nose','python-dateutil','pytz','numpy','astropy',
'sciencedates']
# %%
from setuptools import setup, find_packages

setup(name='dascutils',
      packages=find_packages(),
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/dascutils',
      description='Utilities for UAF Digital All-Sky Camera: reading and plotting',
      version = '0.5',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
      extras_requires={'io':['themisasi'],
                       'plot':['matplotlib','scipy',],},
      install_requires=req,
	  )


