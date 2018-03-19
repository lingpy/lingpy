"""
Setup-Script for LingPy
"""

import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages, Extension
import sys
import os.path
import codecs

# check for specific features
with_c = False
for i, arg in enumerate(sys.argv):
    if arg.lower() == '--with-c':
        del sys.argv[i]
        with_c = True
        break

extra = {}

# setup package name etc as a default
pkgname = 'lingpy'
pkg_dir = {'': '.'}
pkg_location = '.'

if sys.version_info >= (3,):
    extra['use_2to3'] = False

# set up extension modules
if 'install' in sys.argv or 'bdist_egg' in sys.argv or 'develop' in sys.argv:
    if with_c:
        extension_path = ['lingpy', 'algorithm', 'cython']
        extension_prefix = os.path.join(*extension_path)
        extension_modules = [
            Extension(
                os.path.join('.'.join(extension_path + ['calign'])),
                [os.path.join(extension_prefix, 'calign.c')]
            ),
            Extension(
                os.path.join('.'.join(extension_path + ['malign'])),
                [os.path.join(extension_prefix, 'malign.c')]
            ),
            Extension(
                os.path.join('.'.join(extension_path + ['talign'])),
                [os.path.join(extension_prefix, 'talign.c')]
            ),
            Extension(
                os.path.join('.'.join(extension_path + ['cluster'])),
                [os.path.join(extension_prefix, 'cluster.c')]
            ),
            Extension(
                os.path.join('.'.join(extension_path + ['misc'])),
                [os.path.join(extension_prefix, 'misc.c')]
            ),
        ]
    else:
        extension_modules = []
else:
    extension_modules = []

# load requirements
if os.path.isfile('requirements.txt'):
    with open('requirements.txt', 'r') as handle:
        requires = [_.strip() for _ in handle.readlines() if len(_.strip())]
else:
    requires = []

if os.path.isfile('README.md'):
    with codecs.open('README.md', 'r', 'utf-8') as handle:
        long_description=handle.read()
else:
    long_description = ''

# make global name of this version for convenience of modifying it
thisversion = "2.6.3"

setup(
    name=pkgname,
    description="Python library for quantitative tasks in historical linguistics",
    long_description=long_description,
    version=thisversion,
    packages=find_packages(pkg_location, 
        exclude=[ "lingpy._plugins", "_plugins", "*._plugins", "_plugins.*",
            '*._plugins.*', "build", "private", "lingpy.egg-info", "dist",
            "lib"]
        ),
    package_dir=pkg_dir,
    install_requires=requires,
    tests_require=['nose', 'coverage', 'mock'],
    author="Johann-Mattis List and Simon Greenhill and Robert Forkel",
    author_email="info@lingpy.org",
    entry_points={
        'console_scripts' : ['lingpy=lingpy.cli:main'],
    },
    keywords=[
        "historical linguistics",
        "sequence alignment",
        "computational linguistics",
        "dialectology",
        "cognate detection",
    ],
    classifiers=[
        'Natural Language :: English',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Text Processing :: Linguistic',
    ],
    url="http://lingpy.org",
    license="gpl-3.0",
    platforms=["unix", "linux", "windows"],
    ext_modules=extension_modules,
    extras_require={
        "borrowing": ["matplotlib", "scipy"],
        "cldf": ["pycldf"],
        "cluster": ["python-igraph", "scikit-learn"]
    },
    include_package_data=True,
    exclude_package_data={}, 
    package_data={
        '': [
            'data/ipa/sampa.csv',
            'data/swadesh/swadesh.qlc',
            'data/orthography_profiles/*.prf',
            'tests/test_data/*.csv',
            'tests/test_data/*.qlc',
            'tests/test_data/*.msq',
            'tests/test_data/*.msa',
            'tests/test_data/*.tsv',
            'tests/test_data/*.nex',
            'tests/test_data/*.psa', 
            'tests/test_data/*.psq',
            'tests/test_data/cldf/*.csv',
            'tests/test_data/cldf/*.json',
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
            'data/templates/*.html',
            'data/templates/*.js',
            'data/templates/*.css',
            'data/templates/*.tex',
            'data/templates/*.nex',
            'data/swadesh/swadesh.qlc',
        ]
    },
    **extra
    )
