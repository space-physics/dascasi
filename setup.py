#!/usr/bin/env python

from setuptools import setup
import subprocess

try:
    subprocess.call(['conda','install','--file','requirements.txt'])
except Exception as e:
    pass


setup(name='dascutils',
	  description='utilities for the Poker Flat Research Range Digital All Sky Camera, useful for aurora borealis',
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/dascutils',
	  install_requires=['histutils','themisasi',
                        'pathlib2'],
   dependency_links = [
        'https://github.com/scienceopen/histutils/tarball/master#egg=histutils', 
        'https://github.com/scienceopen/themisasi/tarball/master#egg=themisasi'
	],
      packages=['dascutils'],
	  )


