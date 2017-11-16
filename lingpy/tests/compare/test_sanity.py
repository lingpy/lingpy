from __future__ import (
        unicode_literals, print_function, absolute_import, division)

from lingpy.tests.util import test_data
from lingpy.basic.wordlist import Wordlist
from lingpy.compare import sanity as sn

class Tests():
    
    def setUp(self):
        self.wl = Wordlist(test_data('KSL5.qlc'))

    def test__mutual_coverage(self):
        assert len(
                sn._mutual_coverage('Albanian', 'English', self.wl, 'concept')
                ) == 2

    def test__get_concepts(self):
        assert len(
                sn._get_concepts(self.wl, 'concept')) == 6

    def test_mutual_coverage(self):
        assert sn.mutual_coverage(
                self.wl)['French']['Albanian'] == 3

    def test_mutual_coverage_check(self):
        assert not sn.mutual_coverage_check(self.wl, 3)

    def test_mutual_coverage_subset(self):
        a, b = sn.mutual_coverage_subset(
                self.wl, 3, concepts='concept')
        assert a == 3
        assert b[0][0] == 3
        assert b[0][1][0] == 'Albanian'
    
    def test_synonymy(self):
        syns = sn.synonymy(self.wl)
        assert max(syns.values()) == 1
