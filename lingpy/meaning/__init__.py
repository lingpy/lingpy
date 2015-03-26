# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-13 08:26
# modified : 2013-11-13 08:26
"""
Module for the handling of concepts.
"""

__author__="Johann-Mattis List"
__date__="2013-11-13"


from .concepts import ConceptGraph, ConceptComparerBase, \
        ConceptComparerStringMatch, ConceptComparerSpanishStem

from .basvoc import BasVoc

# define a basic basic vocabulary object with the current data we have
concepticon = BasVoc()
