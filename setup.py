# author   : Johann-Mattis List, Peter Bouda, Robert Forkel
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-09 16:28
# modified : 2014-11-05 23:38
"""
Setup-Script for LingPy
"""

__author__="Johann-Mattis List, Peter Bouda, Robert Forkel"
__date__="2014-11-05"


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

# set up extension modules
if 'install' in sys.argv or 'bdist_egg' in sys.argv or 'develop' in sys.argv:
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

requires = [
    'networkx',
    'numpy',
    'six',
    'appdirs'
    ]
    
if sys.version_info < (3, 4):
    requires.append('pathlib')

# make global name of this version for convenience of modifying it
thisversion = "2.4.dev"

setup(
    name=pkgname,
    version=thisversion,
    packages=find_packages(pkg_location),
    package_dir=pkg_dir,
    install_requires=requires,
    tests_require=['regex', 'nltk', 'nose', 'coverage'],
    author="Johann-Mattis List, Steven Moran, Peter Bouda, Johannes Dellert, Taraka Rama, Robert Forkel",
    author_email="mattis.list@uni-marburg.de",
    keywords=[
        "historical linguistics",
        "sequence alignment",
        "computational linguistics"
    ],
    url="http://lingpy.org",
    description="Python library for automatic tasks in historical linguistics",
    license="gpl-3.0",
    platforms=["unix", "linux", "windows"],
    ext_modules=extension_modules,
    extras_require={
        "borrowing": ["matplotlib", "networkx", "scipy"]
    },
    include_package_data=True,
    exclude_package_data={}, 
    package_data={
        '': [
            'data/ipa/sampa.csv',
            'data/orthography_profiles/*.prf',
            'tests/test_data/*.csv',
            'tests/test_data/*.qlc',
            'tests/test_data/*.msq',
            'tests/test_data/*.msa',
            'data/conceptlists/*.tsv',
            'data/conf/*.rc',
            'data/models/*/converter',
            'data/models/*/INFO',
            'data/models/*/matrix',
            'data/models/*/scorer',
            'data/models/dvt/diacritics',
            'data/models/dvt/vowels',
            'data/models/dvt/tones',
            'data/models/dvt_el/diacritics',
            'data/models/dvt_el/vowels',
            'data/models/dvt_el/tones',
            'data/swadesh/swadesh.qlc',
        ]
    },
    **extra)
