# LingPy: A Python Library for Automatic Tasks in Historical Linguistics

This repository contains the Python package `lingpy` which can be used for various tasks in computational historical linguistics.

[![Build Status](https://github.com/lingpy/lingpy/workflows/tests/badge.svg)](https://github.com/lingpy/lingpy/actions?query=workflow%3Atests)
[![codecov.io](http://codecov.io/github/lingpy/lingpy/coverage.svg?branch=master)](http://codecov.io/github/lingpy/lingpy?branch=master)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.597082.svg)](https://doi.org/10.5281/zenodo.597082)
[![PyPI version](https://badge.fury.io/py/lingpy.png)](https://badge.fury.io/py/lingpy)
[![Documentation](https://bade.fury.io/py/lingpy.png)](https://lingpy.github.io)


Authors (Version 2.6.9): Johann-Mattis List and Robert Forkel

Collaborators: 
Christoph Rzymski, Simon J. Greenhill, Steven Moran, Peter Bouda, Johannes Dellert, Taraka Rama, Tiago Tresoldi, Gereon Kaiping, and Frank Nagel.
 
LingPy is a Python library for historical linguistics. It is being developed for Python 2.7 and Python 3.x 
using [a single codebase](https://docs.python.org/3/howto/pyporting.html).

* All source code is available at: [https://github.com/lingpy/lingpy](https://github.com/lingpy/lingpy).
* Documentation can be found at: [http://lingpy.org](http://lingpy.org).
* For a list of papers in which LingPy was applied, see [here](https://github.com/lingpy/lingpy/blob/master/PAPERS.md).

## Quick Installation

For our latest stable version, you can simply use pip or easy_install for installation:
```bash
$ pip install lingpy
```
or 
```bash
$ pip install lingpy
```
Depending on which easy_install or pip version you use, either the Python2 or the Python3 version of LingPy will be installed.

If you want to install the current GitHub version of LingPy on your system, open a terminal and type in the following:
```bash
$ git clone https://github.com/lingpy/lingpy/
$ cd lingpy
$ python setup.py install
```

If the last command above returns you some error regarding user permissions (usually "Errno 13"), you can install
LingPy in your home Python setup:
```
$ python setup.py install --user
```

In order to use the library, start an interactive Python session and import LingPy as follows:
```python
>>> from lingpy import *
```

To install LingPy to hack on it, fork the repository on GitHub, open a terminal and type:
```bash
$ git clone https://github.com/<your-github-user>/lingpy/
$ cd lingpy
$ python setup.py develop
```
This will install LingPy in ["development mode"](http://pythonhosted.org//setuptools/setuptools.html#development-mode),
i.e. you will be able edit the sources in the cloned repository and import the altered code just as the regular Python package.


