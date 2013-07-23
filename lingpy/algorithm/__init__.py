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
from ..settings import rcParams

cmod = {}
# check for c-modules
try:
    from .cython import calign as calign
except ImportError:
    from .cython import _calign as calign
    cmod['calign'] = 1

try:
    from .cython import malign as malign
except ImportError:
    from .cython import _malign as malign
    cmod['malign'] = 1

try:
    from .cython import talign as talign
except:
    from .cython import _talign as talign
    cmod['talign'] = 1

try:
    from .cython import cluster as cluster
except:
    from .cython import _cluster as cluster
    cmod['cluster'] = 1

try:
    from .cython import misc as misc
except:
    from .cython import _misc as misc
    cmod['misc'] = 1

if cmod:
    rcParams['cmodules'] = False
else:
    rcParams['cmodules'] = True
