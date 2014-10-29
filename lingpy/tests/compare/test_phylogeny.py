# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-08-26 10:42
# modified : 2013-08-26 10:42
"""
Test the TreBor borrowing detection algorithm.
"""

__author__="Johann-Mattis List"
__date__="2013-08-26"


import os
from unittest import TestCase

import lingpy
from lingpy.compare.phylogeny import PhyBo
from lingpy.thirdparty.cogent import LoadTree


class TestPhyBo(TestCase):

    def setUp(self):
        self.input = os.path.join(
            os.path.dirname(lingpy.__file__), 'tests', 'test_data', 'phylogeny.qlc')

    def test_get_GLS(self):
        phy = PhyBo(self.input)
        phy.get_GLS()
