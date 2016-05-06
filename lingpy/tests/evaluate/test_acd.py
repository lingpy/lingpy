# coding: utf8
from __future__ import unicode_literals

from lingpy import LexStat
from lingpy.tests.util import test_data, WithTempDir
from lingpy.compare.partial import Partial


class Tests(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.lex = LexStat(test_data('KSL.qlc'))
        self.part = Partial(test_data('partial_cognates.tsv'),
                segments='segments')
        self.part.add_entries('pid1', 'partial_cognate_sets', lambda x: x)
        self.part.add_entries('pid2', 'partialids2', lambda x: [int(y)
            for y in x.split(' ')])


    def test_bcubes(self):
        from lingpy.evaluate.acd import bcubes

        res = bcubes(self.lex, test='cogid', pprint=False)
        self.assertAlmostEquals(res, (1.0, 1.0, 1.0))

        res = bcubes(self.lex, 'cogid', 'cogid', pprint=True,
                per_concept=True)


    def test_partial_bcubes(self):
        from lingpy.evaluate.acd import partial_bcubes
        res = partial_bcubes(self.part, 'pid1', 'pid2', pprint=False)
        assert [round(x, 2) for x in res] == [0.92, 0.98, 0.95]

        res = partial_bcubes(self.part, 'pid1', 'pid2', pprint=True)
        
    def test_pairs(self):
        from lingpy.evaluate.acd import pairs

        res = pairs(self.lex, test='cogid', pprint=False)
        self.assertAlmostEquals(res, (1.0, 1.0, 1.0))

    def test_diff(self):
        from lingpy.evaluate.acd import diff

        res = diff(self.lex, test='cogid', tofile=False, pprint=False)
        self.assertAlmostEquals(res, ((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)))
        self.lex.add_entries('cugid', 'cogid', lambda x: x+1 if x % 2 else x*x)
        d1 = diff(self.lex, gold='cogid', test='cogid', filename='%s' % self.tmp_path('test_acd'), pprint=False)
        d2 = diff(self.lex, gold='cugid', test='cogid', filename='%s' %
                self.tmp_path('test_acd'), pprint=False, tofile=False)
        d3 = diff(self.lex, gold='cugid', test='cogid', filename='%s' %
                self.tmp_path('test_acd'), pprint=False, tofile=True)
        assert d2[0] != 1

def test_npoint_ap():
    from lingpy.evaluate.acd import npoint_ap
    scores = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    scores2 = scores[::-1]
    cognates1 = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    cognates2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    cognates3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    cognates4 = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    np1 = npoint_ap(scores, cognates1)
    np2 = npoint_ap(scores, cognates2)
    np3 = npoint_ap(scores, cognates3)
    np4 = npoint_ap(scores2, cognates1)
    np5 = npoint_ap(scores2, cognates1, reverse=True)
    np6 = npoint_ap(scores, cognates4)
    np7 = npoint_ap(scores2, cognates4)   
    
    assert np2 == 1
    assert np3 == 0
    assert np4 == 0.5
    assert np5 == np1
    assert np6 == 1.0
    assert round(np7, 2) == 0.35 
