# LingPy

[![Build Status](https://travis-ci.org/lingpy/lingpy.svg?branch=master)](https://travis-ci.org/lingpy/lingpy)
[![codecov.io](http://codecov.io/github/lingpy/lingpy/coverage.svg?branch=master)](http://codecov.io/github/lingpy/lingpy?branch=master)


Authors: Johann-Mattis List, Steven Moran, Peter Bouda (Research Unit: [Quantitative Language Comparison](http://www.quanthistling.info/), University of Marburg), and Johannes Dellert (Research Unit: [EVOLAEMP](http://www.sfs.uni-tuebingen.de/~gjaeger/evolaemp/index.html), University of Tübingen), Taraka Rama ([Centre for Language Technology](http://clt.gu.se/), Göteborg), [Robert Forkel](https://github.com/xrotwang) ([Max Planck Institute for the Science of Human History](http://shh.mpg.de/)).

LingPy is a Python Library for Historical Linguistics. It is being developed in Python 3, but we also try to provide basic functionality for Python 2.

* All source code is available at: https://github.com/lingpy/lingpy
* Documentation can be found at: http://lingpy.org 

For a brief overview, see slides from Euroscipy 2012: (http://hprints.org/hprints-00758536)

Other relevant work related to the development of the library and research using it, include:

* List, J. M. 2014. [Sequence comparison in historical linguistics](http://sequencecomparison.github.io). Düsseldorf University Press: Düsseldorf.
* List, J. M. 2012a. Multiple sequence alignment in historical linguistics: a sound class based approach. In Proceedings of Console IXX.
* List, J. M. 2012b. SCA: Phonetic alignment based on sound classes. In New Directions in logic, language, and computation. Slavkovik, M. and D Lassiter (Eds.). Springer.
* Prokic, J. and S. Moran. 2013. Black box approaches to genealogical classification and their shortcomings. Lars Borin and Anju Saxena (eds.) Approaches to Measuring Linguistic Differences. Berlin: Mouton De Gruyter, pages 437–457. 

## Setup for Python 2 and Python 3

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
Note that, although this should now work withouth problems, there is a certain chance that you will need superuser rights the first time you use LingPy. During this first step (which should now already be handled by the setup-script), LingPy compiles some data and stores it as binary pickled objects. Once this has been done, superuser rights are no longer needed.


## Contributing

LingPy development follows the model and workflow outlined here https://gun.io/blog/how-to-github-fork-branch-and-pull-request/


### Setup for Development Version

To install LingPy to hack on it, fork the repository, open a terminal and type:
```bash
$ git clone https://github.com/<your-github-user>/lingpy/
$ cd lingpy
$ python setup.py develop
```
This will just put a symlink to the source code (i.e. to the ``lingpy`` package in the repository clone) in your python site-packages; thus you will be able edit the sources in the git clone and import the altered code just as the regular python package.


### Coding Conventions

In order to keep the code transparent even for multiple contributors, we opt for the following conventions.

#### Importing Modules

When importing external Python modules, try to avoid 

```python
from module import *
```

but instead use

```python
import module
```

In order to avoid heavy typing of long module names, we define the following aliases for thirdparty modules:

```python
import numpy as np
import scipy as sp
import cogent as cg
import networkx as nx
import matplotlip.pyplot as plt
import regex as re # this is useful, since it is then equivalent with the re-module
import mpl_toolkist.basemap as bmp
```

If you are using a specific module offered by numpy or scipy, such as, for example, numpy.linalg, you may import it as:
```python
from numpy import linalg
```
But in order to make it easier for code reviewing to trace thirdparty and internal code, it is recommended to spell out the full names redundantly, unless they occur very often.

Builtin-modules should be imported "as is", without using aliases.

Please make sure to minimize wildcard imports in your code! Wildcard imports are useful in quick-and-dirty scripting, but not in general-purpose libraries. They make it more difficult to trace errors, and can also dramatically increase namespaces. If we load a specific module, we should always know what we want from this module (be it internal or external).

#### Template for Scripts
The first lines of the scripts in which we write the modules should more or less confirm to the following template:

```python
# author   : Max Mustermann
# email    : max.mustermann@uni-muster.com
# created  : 2013-03-14 00:21
# modified : 2013-07-10 11:35
"""
This module provides a basic class for the handling of Mustermanns.

More explicit description with examples can follow here...
"""

__author__="Max Mustermann"
__date__="2013-07-10"

"""
```

- Note that the first not commented line of a scripts should be a short sentence which does not exceed two lines, since sphinx will take this string as a teaser in previews. The following lines can contain a closer description of the respective module.
- Note also that the coding specification which is common in Python2 is not needed in Python3. We add it automatically when converting to Python2.
- The commented line of author, email, etc. is not necessarily required, yet it may turn out usefull to keep a tracker on creation and modification of scripts, and text editors, such as VIM, modify such parts automatically.
- Nevertheless, we should always provide author and date in Python style, just for being able to track down who has done what, and to help people using the library and running into bugs knowing whom to contact.

#### Writing tests

For our tests we use the [nose library](https://nose.readthedocs.org/en/latest/). To run the tests just enter the `lingpy` directory and call `nosetests` or `nosetests.exe` on the command line. *Please do not commit any changes without all tests running without failure or error.*

All our tests are in a directory `tests` within the `lingpy` directory (the latter we will call "source directory" from now on). The `tests` directory mirrors the folder structure of the source directory, i.e. for each directory in the source directory there is a directory in the `tests` directory. For each Python source file in the source directory there is a test file with a prefix "test_". For example, the tests for the module `basic.dictionary`, which has its source in `basic/dictionary.py`, are located in `tests/basic/test_dictionary.py`. Within the test files there is a class defined for each class in the original source files, with a prefix "Test". For example, there is a class `TestDictionary` defined in `test_dictionary.py`, as there is a class `Dictionary` in `dictionary.py`. For each method of a class the test class has a method with the prefix "test_". For example, the method `tokenize()` of the `Dictionary` class is tested with the method `test_tokenize()` of the test class. You may implement a method `setUp()` in the test class to initialize the class you want to test, to load data or do other things that are common for all tests. As an example here is one of the tests of the `Dictionary` class:

```python
import unittest

from lingpy.tests.util import test_data


class TestDictionary(unittest.TestCase):

    def setUp(self):
        self.dictionary = Dictionary(test_data('leach1969-67-161.csv'))

    def test_tokenize(self):
        self.dictionary.tokenize()
        tuples = self.dictionary.get_tuples(['head', 'tokens'])
        assert tuples[0] == ('aa', ['a', 'a'])
```

If all the prefixes are used correctly, then nose will automatically find and execute all the tests.

### Dependencies

Generally, we try to avoid too many dependencies for the creation of LingPy. Without certain modules, however, it 
is difficult to get along without third party libraries. In order to minimize the dependencies,
and to be aware of them ourself, we list them in this section. We distinguish between three different kinds of third
party libraries:

* indispensable libraries, i.e. libraries without which LingPy won't work, 
* recommended libraries, i.e. libraries without which LingPy core functions will work, but which are anyway used quite frequently in a lot of modules, and
* specific libraries, i.e. libraries which serve a specific purpose in a certain marginal module but are not used otherwise

Furthermore, we may conduct temporary forks of libraries that are needed for certain tasks but are not yet available in Python3, 
such as PyCogent (http://pycogent.org/) of which we included some parts in the thirdparty module.

#### Indispensable Libraries

* NumPy: http://numpy.org/
* Regex: http://pypi.python.org/pypi/regex

#### Recommended Libraries

* SciPy: http://scipy.org/
* Networkx: http://networkx.github.com/
* Matplotlib: http://matplotlib.org/

#### Specific Libraries

* Basemap: https://github.com/matplotlib/basemap/

### Creating the Documentation

Currently, documentation is created using the following steps:

* whenever code is added to LingPy, the contributors add documentation inline in their code, following the style also used in other projects such as SciPy, NumPy, and Networkx.
* the general website structure is added around the code, you can find its content by browsing the lingpy/doc/sources/ directory.
* before compiling the code, we create a full reference containing links to all code in LingPy, using the sphinx-apidoc command:

  ```
  $ cd lingpy/ # got to relevenant main folder of lingpy
  $ sphinx-apidoc -o doc/sources/reference/ lingpy
  ```
  
  This creates specific documentation for all the LingPy package, structured with respect to the modules in LingPy.
* run the sphinx build as described on the sphinx website at http://sphinx-doc.org, e.g. on Linux or Mac OSX running

  ```
  $ cd doc
  $ make html
  ```

