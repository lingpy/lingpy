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
import shutil

from three2two import run3to2,run3to3

# check whether a build directory is available
if not os.path.isdir('lingpy_build') and 'install' in sys.argv:
    os.mkdir('lingpy_build')
    os.mkdir(os.path.join('lingpy_build','lingpy'))

# check for install and organize the manifest.in-file
if 'install' in sys.argv:
    m = open('manifest_build.in','r').read()
else:
    m = open('manifest_dist.in','r').read()
f = open('MANIFEST.in','w')
f.write(m)
f.close()


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

    # import lingpy from lingpy_build folder
    if 'install' in sys.argv:
        pkg_location = 'lingpy_build'
        pkgname = 'lingpy'
        pkg_dir = {'':'lingpy_build'}
        run3to3()
        sys.path = ['lingpy_build'] + sys.path
        from lingpy import *
        rc(schema='asjp')
    
else:
    # make a specific directory for lingpy2
    this_version = "2"

    run3to2()
    pkgname = 'lingpy'


    # import lingpy2 to compile the data
    if 'install' in sys.argv:
        # replace manifest path
        pkg_location = 'lingpy_build'
        pkg_dir = {'':'lingpy_build'}
        sys.path = ['lingpy_build'] + sys.path
        from lingpy import *
        rc(schema='asjp')

# set up extension modules
if 'install' in sys.argv or 'bdist_egg' in sys.argv:
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
thisversion = "2.3"
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
if 'install' in sys.argv:
    try:
        print("[i] Removing the build directory.")
        shutil.rmtree('lingpy_build')
        shutil.rmtree('build') 
        print("[i] Done.")
    except:
        print("[i] Failed to remove build directory.")
    print("[i] LingPy was successfully installed on your system.")

