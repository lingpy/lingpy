"""
Test the TreBor borrowing detection algorithm.
"""
from unittest import TestCase

from lingpy.compare.phylogeny import PhyBo
from lingpy.tests.util import test_data


class TestPhyBo(TestCase):
    def setUp(self):
        self.phy = PhyBo(test_data('sagart.qlc'))
