# author   : Johann-Mattis List, Peter Bouda
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-09 16:28
# modified : 2013-10-18 11:35
"""
Setup-Script for LingPy
"""

__author__="Johann-Mattis List,Peter Bouda"
__date__="2013-10-18"


import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages,Extension
import sys
import os
import os.path
import shutil

from three2two import run3to2,run3to3

# check whether a build directory is available
if not os.path.isdir('lingpy_build'):
    os.mkdir('lingpy_build')
    os.mkdir('lingpy_build/lingpy')

# check for specific features
with_c = False
for i,arg in enumerate(sys.argv):
    if arg.lower() == '--with-c':
        del sys.argv[i]
        with_c = True
        break

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = False
    this_version = "3"
    pkg_location = 'lingpy_build'
    pkgname = 'lingpy'
    pkg_dir = {'':'lingpy_build'}
    run3to3()

    # import lingpy from lingpy_build folder
    sys.path = ['lingpy_build/'] + sys.path
    from lingpy import *
    rc(schema='asjp')
    
else:
    # make a specific directory for lingpy2
    this_version = "2"

    run3to2()
    pkgname = 'lingpy'
    # replace manifest path
    pkg_location = 'lingpy_build'
    pkg_dir = {'':'lingpy_build'}

    # import lingpy2 to compile the data
    sys.path = ['lingpy_build/'] + sys.path
    from lingpy import *
    rc(schema='asjp')

# set up extension modules
if 'install' in sys.argv or 'bdist_egg' in sys.argv:
    if this_version == "3":

        if with_c: #"--with-c" in sys.argv or '--with-C' in sys.argv in sys.argv:
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
thisversion = "2.1"
setup(
        name = pkgname,
        version = thisversion,
        packages = find_packages(pkg_location),
        package_dir = pkg_dir,
        install_requires = ['numpy','regex'],
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

# remove the build directory in order to prevent that it leads to confusion
# when installing lingpy in both the py2 and the py3 version
print("[i] Removing the build directory.")
shutil.rmtree('lingpy_build/')
print("[i] Done.")
print("[i] LingPy was successfully installed on your system.")
