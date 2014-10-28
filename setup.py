# author   : Johann-Mattis List, Peter Bouda
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-09 16:28
# modified : 2013-11-07 13:01
"""
Setup-Script for LingPy
"""

__author__="Johann-Mattis List,Peter Bouda"
__date__="2013-11-07"


import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages,Extension
import sys
import os
import os.path


# check for specific features
with_c = False
for i,arg in enumerate(sys.argv):
    if arg.lower() == '--with-c':
        del sys.argv[i]
        with_c = True
        break

extra = {}

# setup package name etc as a default
pkgname = 'lingpy'
pkg_dir = {'':'.'}
pkg_location = '.'

if sys.version_info >= (3,):
    extra['use_2to3'] = False
    this_version = "3"
else:
    # make a specific directory for lingpy2
    this_version = "2"

from lingpy import *
rc(schema='asjp')

# set up extension modules
if 'install' in sys.argv or 'bdist_egg' in sys.argv or 'develop' in sys.argv:
    if this_version == "3":

        if with_c:
            extension_path = ['lingpy','algorithm','cython']
            extension_prefix = os.path.join(*extension_path)
            extension_modules = [
                        Extension(
                            os.path.join('.'.join(extension_path),'calign'),
                            [os.path.join(extension_prefix,'calign.c')]
                            ),
                        Extension(
                             os.path.join('.'.join(extension_path),'malign'),
                            [os.path.join(extension_prefix,'malign.c')]
                            ),
                        Extension(
                            os.path.join('.'.join(extension_path),'talign'),
                            [os.path.join(extension_prefix,'talign.c')]
                            ),
                        Extension(
                            os.path.join('.'.join(extension_path),'cluster'),
                            [os.path.join(extension_prefix,'cluster.c')]
                            ),
                        Extension(
                            os.path.join('.'.join(extension_path),'misc'),
                            [os.path.join(extension_prefix,'misc.c')]
                            ),
                        ]
        else:
            extension_modules = []

    else:

        extension_modules = []


else:
    extension_modules = []

# make global name of this version
thisversion = "2.4.dev"
setup(
        name = pkgname,
        version = thisversion,
        packages = find_packages(pkg_location),
        package_dir = pkg_dir,
        install_requires = [
		'numpy', 
		'six',
		'networkx',
		'regex',
		#'matplotlib',
		#'scipy',
	],
        tests_require=['nltk', 'nose', 'coverage'],
        author = "Johann-Mattis List, Steven Moran, Peter Bouda, Johannes Dellert",
        author_email = "mattis.list@uni-marburg.de,steven.moran@lmu.de",
        keywords = [
            "historical linguistics", 
            "sequence alignment",
            "computational linguistics"
            ],
        url = "http://lingpy.org",
        description = "Python library for automatic tasks in historical linguistics",
        license = "gpl-3.0",
        platforms = ["unix","linux","windows"],
        ext_modules=extension_modules,
        extras_require = {
            "borrowing" : ["matplotlib","networkx","scipy"]
            },
        include_package_data = True,
        exclude_package_data = {}, #{'':["*.bin"]},
        **extra
        )

