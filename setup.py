# LingPy 
#
# Copyright (C) 2012 Johann-Mattis List 
# Author: Johann-Mattis List, Steven Moran
# Author Email: <mattis.list@uni-marburg.de> URL: <http://lingulist.de/lingpy> For
# license information, see <gpl-3.0.txt>

import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages,Extension
import sys

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = False
    #extra['convert_2to3_doctests'] = ['src/your/module/README.txt']
    #extra['use_2to3_fixers'] = ['your.fixers']

setup(
        name = "lingpy",
        version = "2.0.dev",
        packages = find_packages(),
        include_package_data = True,
        install_requires = ['numpy','networkx'],
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
        platforms = ["unix","linux"],
        ext_modules=[
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
            ],
        **extra
        )


