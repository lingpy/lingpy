# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-07 16:56
# modified : 2013-03-07 16:56
"""
Package provides basic modules for alignment analyses.
"""

__author__="Johann-Mattis List"
__date__="2013-03-07"


from .multiple import Multiple, mult_align
from .pairwise import Pairwise, pw_align, nw_align, sw_align, we_align, \
        structalign, turchin, edit_dist
from .sca import MSA, PSA, Alignments, SCA
