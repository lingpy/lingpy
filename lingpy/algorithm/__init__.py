"""
Package for specific algorithms and time-intensive routines.
"""
# flake8: noqa
from lingpy.settings import rcParams
from lingpy.algorithm.clustering import *
from lingpy.algorithm._tree import _TreeDist as TreeDist

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
    from .cython import misc as misc
except:
    from .cython import _misc as misc
    cmod['misc'] = 1

rcParams['cmodules'] = not bool(cmod)

# define squareform for global lingpy-applications
squareform = misc.squareform
