# -*- coding: UTF-8 -*-
#! /usr/bin/python

from pathlib import Path
from setuptools import setup, find_packages

# ...
# Read library version into '__version__' variable
path = Path(__file__).parent / 'PARPY6P' / 'version.py'
exec(path.read_text())
# ...

NAME    = 'PARPY6P'
VERSION = __version__
AUTHOR  = 'Imad Kissami'
EMAIL   = 'imad.kissami@gmail.com'
URL     = 'https://github.com/imadki/PARPY6P'
DESCR   = 'Python extension language using accelerators.'
KEYWORDS = ['math']
LICENSE = "LICENSE"

setup_args = dict(
    name                 = NAME,
    version              = VERSION,
    description          = DESCR,
    #long_description     = open('README.rst').read(),
    author               = AUTHOR,
    author_email         = EMAIL,
    license              = LICENSE,
    keywords             = KEYWORDS,
#    url                  = URL,
)

# ...
packages = find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"])
# ...

# Dependencies
install_requires = [
    'numpy',
    'meshio',
    'mpi4py',
    'metis',
    #'timeit',
    ]

def setup_package():
    setup(packages=packages, \
          include_package_data=True, \
          install_requires=install_requires, \
         **setup_args)

if __name__ == "__main__":
    setup_package()
