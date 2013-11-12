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
from lingpy.compare.phylogeny import PhyBo

class TestPhyBo:

    def setup(self):

        self.phy = PhyBo('.../test_data/sagart.qlc')

