#!/usr/bin/env python3

from setuptools import setup
import subprocess

try:
    subprocess.run(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    print('you will need to install packages in requirements.txt  {}'.format(e))


with open('README.rst','r') as f:
	long_description = f.read()

setup(name='dascutils',
      version='0.1',
	  description='utilities for the Poker Flat Research Range Digital All Sky Camera, useful for aurora borealis',
	  long_description=long_description,
	  author='Michael Hirsch',
	  url='https://github.com/scienceopen/dascutils',
	  install_requires=['histutils','themisasi'],
   dependency_links = [
        'https://github.com/scienceopen/histutils/tarball/master#egg=histutils', 
        'https://github.com/scienceopen/themisasi/tarball/master#egg=themisasi'],
      packages=['dascutils'],
	  )


