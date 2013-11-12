# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-12 12:46
# modified : 2013-11-12 12:46
"""
Test lexstat module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-12"

import os
from lingpy import LexStat
from lingpy.settings import rcParams

class TestLexStat:

    def setup(self):

        self.lex = LexStat(
                os.path.join(
                    rcParams['_path'],
                    'tests',
                    'test_data',
                    'KSL.qlc'
                    )
                )

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

