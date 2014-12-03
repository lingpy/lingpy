# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-12 12:46
# modified : 2013-11-12 12:46
"""
Test lexstat module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-12"

import unittest
from lingpy import LexStat
from lingpy.tests.util import test_data


class TestLexStat(unittest.TestCase):

    def setUp(self):
        self.lex = LexStat(test_data('KSL.qlc'))

    def test_get_scorer(self):
        self.lex.get_scorer()
        assert hasattr(self.lex,"cscorer") == True

    def test_cluster(self):
        self.lex.get_scorer()
        self.lex.cluster(method="lexstat", threshold=0.7)
        self.lex.cluster(method="edit-dist", threshold=0.7)
        self.lex.cluster(method="turchin", threshold=0.7)

        assert ('scaid' in self.lex.header and 'lexstatid' in self.lex.header \
                and 'editid' in self.lex.header and 'turchinid' in \
                self.lex.header) == True
    
    def test_align_pairs(self):
        self.lex.align_pairs('English', 'German', method='sca')
