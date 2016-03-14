# Coding Convetions for LingPy

First and foremost: Follow [PEP8](https://www.python.org/dev/peps/pep-0008/).


## General module layout

LingPy python modules should conform to the following template:

```python
# coding: utf8
"""
This module provides a basic class for the handling of Mustermanns.

More explicit description with examples can follow here...
"""
from __future__ import unicode_literals, print_function, division
```

The docstring of a module should consist of a short sentence which does 
not exceed two lines, since sphinx will take this string as a teaser in previews. The 
following lines can contain a closer description of the respective module.


## Importing Modules

The following common aliases for thirdparty modules can be used in `import`s:

```python
import numpy as np
import scipy as sp
import cogent as cg
import networkx as nx
import matplotlip.pyplot as plt
import regex as re  # this is useful, since it is then equivalent with the re-module
import mpl_toolkist.basemap as bmp
```

If you are using a specific module offered by numpy or scipy, such as, for example, numpy.linalg, you may import it as:
```python
from numpy import linalg
```
But in order to make it easier for code reviewing to trace thirdparty and internal code, 
it is recommended to spell out the full names redundantly, unless they occur very often.

Builtin-modules should be imported "as is", without using aliases.

Please make sure to use no wildcard imports in your code! Wildcard imports are useful in 
quick-and-dirty scripting, but not in general-purpose libraries. They make it more difficult 
to trace errors, and can also dramatically increase namespaces. If we load a specific module, 
we should always know what we want from this module (be it internal or external).


## Tests

All our tests are in a directory `tests` within the `lingpy` package (the latter we will 
call "source directory" from now on). The `tests` directory mirrors the folder structure 
of the source directory, i.e. for each directory in the source directory there is a 
directory in the `tests` directory. For each Python module in the source directory there 
is a test file with a prefix "test_". For example, the tests for the module 
`basic.dictionary`, which has its source in `basic/dictionary.py`, are located in 
`tests/basic/test_dictionary.py`. 

Within the test files there is a class defined for each 
class in the original source files, with a prefix "Test". For example, there is a class 
`TestDictionary` defined in `test_dictionary.py`, as there is a class `Dictionary` in 
`dictionary.py`. For each method of a class the test class has a method with the prefix 
"test_". For example, the method `tokenize()` of the `Dictionary` class is tested with the 
method `test_tokenize()` of the test class. 

All tests should be implemented as subclasses of [`unittest.testCase`](https://docs.python.org/3/library/unittest.html).
Thus, you may implement a method `setUp()` in the 
test class to initialize the class you want to test, to load data or do other things that 
are common for all tests. As an example here is one of the tests of the `Dictionary` class:

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

If functionality excuted in a test creates files in the filesystem, the test class must 
inherit from `lingpy.tests.util.WithTmpDir` and care must be taken to only create files
in the temporary directory which is provided (and disposed off) by this test class.


## Logging

LingPy uses [python's `logging` module](https://docs.python.org/3/library/logging.html) 
for logging. To make sure an appropriately configured logger is used, library code must
use the functionality in `lingpy.log` to get a logger instance or emit log messages.

Use of the `print` function must always be guarded by a check for the log level:
```python
import logging
from lingpy import log

...

    if log.get_level() <= logging.INFO:
        print('...')
```


## Dependencies

Generally, we try to avoid too many dependencies for the creation of LingPy. Without certain modules, however, it 
is difficult to get along without third party libraries. In order to minimize the dependencies,
and to be aware of them ourself, we list them in this section. We distinguish between three different kinds of third
party libraries:

* indispensable libraries, i.e. libraries without which LingPy won't work, 
* recommended libraries, i.e. libraries without which LingPy core functions will work, but which are anyway used quite frequently in a lot of modules, and
* specific libraries, i.e. libraries which serve a specific purpose in a certain marginal module but are not used otherwise

Furthermore, we may conduct temporary forks of libraries that are needed for certain tasks but are not yet available in Python3, 
such as PyCogent (http://pycogent.org/) of which we included some parts in the thirdparty module.

### Indispensable Libraries

* NumPy: http://numpy.org/

### Recommended Libraries

* SciPy: http://scipy.org/
* Networkx: http://networkx.github.com/
* Matplotlib: http://matplotlib.org/

### Specific Libraries

* Basemap: https://github.com/matplotlib/basemap/
