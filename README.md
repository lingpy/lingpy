# QLC-LingPy

Authors: Johann-Mattis List, Steven Moran, Peter Bouda (Research Unit: [Quantitative Language Comparison](http://www.quanthistling.info/), University of Marburg), and Johannes Dellert (Research Unit: [EVOLAEMP](http://www.sfs.uni-tuebingen.de/~gjaeger/evolaemp/index.html), University of Tübingen).

QLC-LingPy is a Python Library for Historical Linguistics. It is being developed in Python 3, but we also try to provide basic functionality for Python 2.

Available at: (https://github.com/lingpy/lingpy)

For a brief overview, see slides from Euroscipy 2012: (http://hprints.org/hprints-00758536)

Other relevant work related to the development of the library and research using it, include:

* List, J. M. Forthcoming. Sequence comparison in historical linguistics. PhD dissertation, University of Düsseldorf.
* List, J. M. 2012a. Multiple sequence alignment in historical linguistics: a sound class based approach. In Proceedings of Console IXX.
* List, J. M. 2012b. SCA: Phonetic alignment based on sound classes. In New Directions in logic, language, and computation. Slavkovik, M. and D Lassiter (Eds.). Springer.
* Prokic, J. and S. Moran. 2012. Black box approaches to genealogical classication and their shortcomings. To appear in Comparing Approaches to Measuring Linguistic Dierences, Lars Borin and Anju Saxena (eds). Mouton De Gruyter. 


## Setup

To use the library in its pre-setup.py stage, git clone the library and put the library's "lingpy/lingpy" folder in your path (or Python path).

Alternatively, you can make a symlink in your Python 3 site-packages folder called "lingpy" and link it to "lingpy/lingpy". For example:

1. Start the Python interpreter (make sure you are using Python 3)

```bash
$ python
```

2. At the prompt, locate the site-packages folder:

```python
>>> import sys
>>> print(sys.path)
['', '/opt/local/Library/Frameworks/Python.framework/Versions/3.2/lib/python3.2', '/opt/local/Library/Frameworks/Python.framework/Versions/3.2/lib/python3.2/site-packages']
```

3. Create the symlink (you may need sudo):

```bash
ln -s /path/to/lingpy/lingpy /path/to/site-packages/symlink
```

4. Test in interpreter:

```bash
$ python
```

```python
>>> import lingpy
```
## Coding Conventions

In order to keep the code transparent even for multiple contributors, we opt for the following conventions.

### Importing Modules

When importing Python modules, try to avoid 

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
```

Builtin-modules should be imported "as is", without using aliases.

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
* Regex: http://pypi.python.org/pypi/regex

### Recommended Libraries

* SciPy: http://scipy.org/
* Networkx: http://networkx.github.com/
* Matplotlib: http://matplotlib.org/

### Specific Libraries

* Basemap: https://github.com/matplotlib/basemap/






