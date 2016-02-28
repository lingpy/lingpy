# LingPy

[![Build Status](https://travis-ci.org/lingpy/lingpy.svg?branch=master)](https://travis-ci.org/lingpy/lingpy)
[![codecov.io](http://codecov.io/github/lingpy/lingpy/coverage.svg?branch=master)](http://codecov.io/github/lingpy/lingpy?branch=master)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.16093.svg)](http://dx.doi.org/10.5281/zenodo.16093)
[![Latest Stable Version](https://poser.pugx.org/phpunit/phpunit/version)]https://github.com/lingpy/lingpy/releases/tag/v2.3)
[![Latest Unstable Version](https://poser.pugx.org/phpunit/phpunit/v/unstable)](https://github.com/lingpy/lingpy/releases/tag/v2.4.1-alpha)


Authors: [Johann-Mattis List](https://github.com/linguist) ([CRLAO, Paris](http://crlao.ehess.fr/)) and [Robert Forkel](https://github.com/xrotwang) ([Max Planck Institute for the Science of Human History](http://shh.mpg.de/))

Collaborators: Steven Moran ([Universität Zürich](http://www.linguistik.uzh.ch/about/mitglieder/moran.html), [Peter Bouda](http://www.peterbouda.eu/), Johannes Dellert ([University of Tübingen](http://www.sfs.uni-tuebingen.de/~gjaeger/evolaemp/index.html)), Taraka Rama ([Centre for Language Technology](http://clt.gu.se/), Göteborg), Simon Greenhill ([Australian National University, Canberra](https://researchers.anu.edu.au/researchers/greenhill-s).

LingPy is a Python Library for Historical Linguistics. It is being developed in Python 3, but we also provide basic functionality for Python 2.

* All source code is available at: [https://github.com/lingpy/lingpy](https://github.com/lingpy/lingpy).
* Documentation can be found at: [http://lingpy.org](http://lingpy.org).
* For a list of papers in which LingPy was applied, see [here](https://github.com/lingpy/lingpy/blob/master/PAPERS.md).

## Quick Installation

If you want to regularly install LingPy on your system, open a terminal and type in the following:
```bash
$ git clone https://github.com/lingpy/lingpy/
$ cd lingpy
$ python setup.py install
```

In order to use the library, open Python2 or Python3 in your terminal and import LingPy as follows:
```python
>>> from lingpy import *
```
To install LingPy to hack on it, fork the repository, open a terminal and type:
```bash
$ git clone https://github.com/<your-github-user>/lingpy/
$ cd lingpy
$ python setup.py develop
```
This will just put a symlink to the source code (i.e. to the ``lingpy`` package in the repository clone) in your python site-packages; thus you will be able edit the sources in the git clone and import the altered code just as the regular python package.


## Contributing

LingPy development follows the model and workflow outlined here https://gun.io/blog/how-to-github-fork-branch-and-pull-request/. In order to keep the code transparent even for multiple contributors, we have set up a list of coding conventions. See [here](https://github.com/lingpy/lingpy/blob/master/CONVENTIONS.md) for details.


