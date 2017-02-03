#!/usr/bin/env python
from setuptools import setup

try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)


setup(name='dascutils',
	  install_requires=['histutils','themisasi',],
     dependency_links = [
        'https://github.com/scienceopen/histutils/tarball/master#egg=histutils', 
        'https://github.com/scienceopen/themisasi/tarball/master#egg=themisasi'
	],
      packages=['dascutils'],
	  )


