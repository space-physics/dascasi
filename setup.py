#!/usr/bin/env python
install_requires = ['python-dateutil','numpy>=1.13','astropy','xarray']
tests_require=['pytest','nose','coveralls']
# %%
from setuptools import setup, find_packages

setup(name='dascutils',
      packages=find_packages(),
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/dascutils',
      description='Utilities for UAF Digital All-Sky Camera: reading and plotting',
      long_description=open('README.rst').read(),
      version = '1.2.1',
      classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3.6',
      'Programming Language :: Python :: 3.7',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      extras_require={'io':['themisasi'],
                       'plot':['matplotlib','scipy',],
                       'tests':tests_require},
      install_requires=install_requires,
      tests_require=tests_require,
      python_requires='>=3.6',
      scripts=['DownloadDASC.py','PlotDASC.py']
	  )


