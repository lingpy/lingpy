# author   : Johann-Mattis List, Steven Moran
# email    : mattis.list@gmail.com
# created  : 2013-03-04 14:05
# modified : 2014-08-15 13:39
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

__author__="Johann-Mattis List, Steven Moran, Peter Bouda, Johannes Dellert"
__date__="2014-08-15"

# import settings
from .settings import rc

# general imports
from .basic import Wordlist, Dictionary, Tree, Workflow

# import converts
# from .convert import *

# we don't import align for the moment for safety reasons...
from .align import Pairwise, Multiple, SCA, MSA, PSA, Alignments, edit_dist, \
        pw_align, nw_align, sw_align, we_align, structalign, turchin, \
        mult_align

# load the sound-class models
from .data import Model

# import reading routine
from .read import csv2list, csv2dict, star2qlc

# import sequence routines
from .sequence import ipa2tokens, tokens2class, prosodic_string, \
        prosodic_weights, class2tokens, pid, sampa2uni

# import thirdparty modules
#from .thirdparty import LoadTree as Tree

# import lexstat
from .compare.lexstat import LexStat

# import algorithm-stuff
from .algorithm.clustering import upgma, neighbor, flat_upgma, flat_cluster, \
        fuzzy, link_clustering, mcl

