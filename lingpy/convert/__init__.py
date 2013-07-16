# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-06-11 22:50
# modified : 2013-06-11 22:50
"""
Package provides different methods for file conversion.
"""

__author__="Johann-Mattis List"
__date__="2013-06-11"


from .nexus import pap2nex
from .csv import pap2csv,wl2csv
from .phylip import matrix2dst
from .misc import *
from .newick import *
from .utils import *
