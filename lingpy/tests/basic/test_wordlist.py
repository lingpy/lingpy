# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-12 10:27
# modified : 2013-11-12 10:27
"""
Test wordlist module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-12"

import os
from lingpy import Wordlist
from lingpy.settings import rcParams

class TestWordlist:

    def setup(self):
        self.wordlist = Wordlist(
                os.path.join(
                    rcParams['_path'],
                    'tests',
                    'test_data',
                    'KSL.qlc'
                    )
                )

    def test___len__(self):

        assert len(self.wordlist) == 1400

    def test_calculate(self):

        self.wordlist.calculate('tree')

        assert sorted(self.wordlist.tree.taxa) == sorted(self.wordlist.cols)

        


