# author   : Johann-Mattis List, Steven Moran
# email    : mattis.list@gmail.com
# created  : 2013-03-04 14:05
# modified : 2013-04-04 22:48
"""
LingPy package for quantitative tasks in historical linguistics.

Documentation is available in the docstrings. Online documentation is available
at http://lingpy.org

Subpackages
-----------
algorithm  --- Basic Algorithms for Sequence Comparison
align      --- Specific Algorithms Alignment Analyses
basic      --- Basic Classes for Language Comparison
check      --- Classes for Exceptions, Warnings, and Check
compare    --- Basic Modules for Language Comparison
convert    --- Functions for Format Conversion
data       --- Data Handling
evaluate   --- Basic Classes and Functions for Algorithm Evaluation
read       --- Basic Functions for Data Input
sequence   --- Basic Functions for Sequence Modeling
thirdparty --- Temporary Forks of Third-Party-Modules

"""

__author__="Johann-Mattis List, Steven Moran"
__date__="2013-04-04"

# general imports
from .basic import *

# import converts
from .convert import *

# we don't import align for the moment for safety reasons...
from .align import *

# load the sound-class models
from .data import *

# import reading routines
from .read import *

# import sequence routines
from .sequence import *

# import thirdparty modules
from .thirdparty import *

# import lexstat
from .compare.lexstat import *

# import algorithm-stuff
from .algorithm import *
