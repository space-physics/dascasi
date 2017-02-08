#!/usr/bin/env python
from setuptools import setup

try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)
    import pip
    pip.main(['install','-r','requirements.txt'])

setup(name='dascutils',
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scienceopen/dascutils',
      description='Utilities for UAF Digital All-Sky Camera: reading and plotting',
      version = '0.5',
      classifiers=[
      'Programming Language :: Python :: 3.6',
      ],
	  install_requires=['sciencedates',
	                    'themisasi',],
     dependency_links = [
        'https://github.com/scienceopen/themisasi/tarball/master#egg=themisasi'
	],
      packages=['dascutils'],
	  )


