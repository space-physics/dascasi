#!/usr/bin/env python
install_requires = ['python-dateutil','pytz','numpy','astropy',
'sciencedates']
tests_require=['nose','coveralls']
# %%
from setuptools import setup, find_packages

setup(name='dascutils',
      packages=find_packages(),
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/dascutils',
      description='Utilities for UAF Digital All-Sky Camera: reading and plotting',
      long_description=open('README.rst').read(),
      version = '0.5.1',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
      extras_requires={'io':['themisasi'],
                       'plot':['matplotlib','scipy',],},
      install_requires=install_requires,
      tests_require=tests_require,
      python_requires='>=3.6',
	  )


