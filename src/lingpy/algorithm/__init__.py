"""
Package for specific algorithms and time-intensive routines.
"""
# flake8: noqa
from lingpy.settings import rcParams
from lingpy.algorithm.clustering import *
from lingpy.algorithm._tree import _TreeDist as TreeDist

from .cython import _calign as calign
from .cython import _malign as malign
from .cython import _talign as talign
from .cython import _misc as misc

rcParams['cmodules'] = False

# define squareform for global lingpy-applications
squareform = misc.squareform
