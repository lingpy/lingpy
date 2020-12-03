"""
Package provides basic modules for alignment analyses.
"""
# flake8: noqa
from .multiple import Multiple, mult_align
from .pairwise import (
    Pairwise, pw_align, nw_align, sw_align, we_align, structalign, turchin, edit_dist,
)
from .sca import MSA, PSA, Alignments, SCA
