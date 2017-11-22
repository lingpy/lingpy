# LingPy

[![Build Status](https://travis-ci.org/lingpy/lingpy.svg?branch=master)](https://travis-ci.org/lingpy/lingpy)
[![codecov.io](http://codecov.io/github/lingpy/lingpy/coverage.svg?branch=master)](http://codecov.io/github/lingpy/lingpy?branch=master)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.597082.svg)](https://doi.org/10.5281/zenodo.597082)
[![PyPI version](https://badge.fury.io/py/lingpy.png)](https://badge.fury.io/py/lingpy)
[![Documentation](https://bade.fury.io/py/lingpy.png)](https://lingpy.github.io)


Authors: [Johann-Mattis List](https://github.com/linguist) ([Max Planck Institute for the Science of Human History](http://shh.mpg.de/)), [Simon Greenhill](https://github.com/simongreenhill) ([Max Planck Institute for the Science of Human History](http://shh.mpg.de/)), and [Robert Forkel](https://github.com/xrotwang) ([Max Planck Institute for the Science of Human History](http://shh.mpg.de/))

Collaborators: 
Steven Moran ([Universität Zürich](http://www.linguistik.uzh.ch/about/mitglieder/moran.html), [Peter Bouda](http://www.peterbouda.eu/), Johannes Dellert ([University of Tübingen](http://www.sfs.uni-tuebingen.de/~gjaeger/evolaemp/index.html)), Taraka Rama ([Centre for Language Technology](http://clt.gu.se/), Göteborg), and Tiago Tresoldi ([UFRGS](http://www.ufrgs.br/english/home))
 
LingPy is a Python Library for Historical Linguistics. It is being developed for Python 2.7 and Python 3.x 
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
$ easy_install lingpy
```
Depending on which easy_install or pip version you use, either the Python2 or the Python3 version of LingPy will be installed.

If you want to install LingPy the current GitHub version on your system, open a terminal and type in the following:
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

In order to use the library, start an interactive python session and import LingPy as follows:
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
i.e. you will be able edit the sources in the cloned repository and import the altered code just as the regular python package.


