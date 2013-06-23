# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-12 08:05
# modified : 2013-04-15 09:30
"""
Package for specific algorithms and time-intensive routines.
"""

__author__="Johann-Mattis List"
__date__="2013-04-15"

from .distance import *

# check for c-modules
try:
    from .cython import calign as calign
    from .cython import malign as malign
    from .cython import talign as talign
    from .cython import cluster as cluster
    from .cython import misc as misc
except:
    print(
        "[i] Import of C-modules failed, using pure Python implementation instead."
        )
    from .cython import _calign as calign
    from .cython import _malign as malign
    from .cython import _talign as talign
    from .cython import _cluster as cluster
    from .cython import _misc as misc


