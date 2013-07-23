# LingPy 
#
# Copyright (C) 2012 Johann-Mattis List 
# Author: Johann-Mattis List, Steven Moran
# Author Email: <mattis.list@uni-marburg.de> URL: <http://lingulist.de/lingpy> 
# For license information, see <gpl-3.0.txt>

import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages,Extension
import sys
import os
import os.path

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = False
    this_version = "3"
    pkgname = 'lingpy'
else:
    # make a specific directory for lingpy2
    this_version = "2"
    if not os.path.isdir('lingpy2'):
        os.mkdir('lingpy2')
    from three2two import run3to2
    run3to2()

    pkgname = 'lingpy2'
    # replace manifest path

if this_version == '2':
    f = open('MANIFEST.in').read()
    if not 'lingpy2' in f:
        out = open('MANIFEST.in','w')
        out.write(f.replace('lingpy','lingpy2'))
        out.close()
else:
    f = open('MANIFEST.in').read()
    if 'lingpy2' in f:
        out = open('MANIFEST.in','w')
        out.write(f.replace('lingpy2','lingpy'))
        out.close()


# set up extension modules
if 'install' in sys.argv or 'bdist_egg' in sys.argv:
    if this_version == "3":
        answer = input("[i] Do you want to install with C-modules (requires Cython)? (Y/N) ")

        if answer.lower() in ['y','j','yes']:
            extension_modules = [
                        Extension(
                            'lingpy.algorithm.cython/calign',
                            ['lingpy/algorithm/cython/calign.c'],
                            ),
                        Extension(
                            'lingpy.algorithm.cython/malign',
                            ['lingpy/algorithm/cython/malign.c'],
                            ),
                        Extension(
                            'lingpy.algorithm.cython/talign',
                            ['lingpy/algorithm/cython/talign.c'],
                            ),
                        Extension(
                            'lingpy.algorithm.cython/cluster',
                            ['lingpy/algorithm/cython/cluster.c'],
                            ),
                        Extension(
                            'lingpy.algorithm.cython/misc',
                            ['lingpy/algorithm/cython/misc.c'],
                            ),
                        ]
        else:
            extension_modules = []

    else:

        extension_modules = []


else:
    extension_modules = []

# make global name of this version
thisversion = "2.0"
setup(
        name = pkgname,
        version = thisversion,
        packages = find_packages(),
        include_package_data = True,
        install_requires = ['numpy','networkx','regex'],
        author = "Johann-Mattis List, Steven Moran",
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
        **extra
        )

if this_version == '2':
    from lingpy2 import *
else:
    from lingpy import *
