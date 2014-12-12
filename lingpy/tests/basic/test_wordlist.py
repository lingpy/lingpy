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

from six import text_type

from lingpy import Wordlist
from lingpy.settings import rcParams
from lingpy.tests.util import test_data, WithTempDir


class TestWordlist(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.wordlist = Wordlist(test_data('KSL.qlc'))

    def test___len__(self):
        assert len(self.wordlist) == 1400

    def test_calculate(self):
        self.wordlist.calculate('dst')

        assert hasattr(self.wordlist,'distances')
        assert sum([self.wordlist.distances[x][x] for x in
            range(self.wordlist.width)]) == 0

        self.wordlist.calculate('tree')

        assert sorted(self.wordlist.tree.taxa) == sorted(self.wordlist.cols)

        self.wordlist.calculate('groups')

        assert hasattr(self.wordlist,'groups')
        assert type(self.wordlist.groups) == dict

    def test_get_list(self):
        gerL = self.wordlist.get_list(col='German', entry='ipa', flat=True)
        gerD = self.wordlist.get_dict(col='German',entry='ipa')

        assert sorted(gerL) == sorted([v[0] for v in gerD.values()])

    def test_get_dict(self):
        gerD = self.wordlist.get_dict(col='German')

        assert sorted(gerD.keys()) == sorted(self.wordlist.rows)

    def test_renumber(self):
        self.wordlist.renumber('cogid','dummy')

        ger1 = self.wordlist.get_list(col='German', entry='cogid', flat=True)
        ger2 = self.wordlist.get_list(col='German', entry='dummy', flat=True)

        assert len(set(ger1)) == len(set(ger2))
        assert sum([1 for x in ger2 if type(x) == int]) == len(ger2)

    def test_get_entries(self):
        ger = self.wordlist.get_entries('cogid')

        assert len(ger) == self.wordlist.height
        assert len(ger[0]) == self.wordlist.width

    def get_etymdict(self):
        etd1 = self.wordlist.get_etymdict(ref='cogid', entry='ipa',
                loans=False)
        etd2 = self.wordlist.get_etymdict(ref='cogid', entry='ipa',
                loans=True)

        assert len(etd1) > len(etd2) and len(set([abs(x) for x in etd1])) == \
                len(etd2)
        assert len([x for x in etd2 if x < 0]) == 0

        # make "fuzzy" cognate sets
        self.wordlist.add_entries(
                'fuzzyid',
                'cogid',
                lambda x: [x]
                )

        etd3 = self.wordlist.get_etymdict(ref='fuzzyid', entry='ipa',
                loans=False, fuzzy=True)
        etd4 = self.wordlist.get_etymdict(ref='fuzzyid', entry='ipa',
                loans=True, fuzzy=True)
        for key in etd1:
            assert etd1[key] == etd3[key]
        for key in etd2:
            assert etd2[key] == etd4[key]

    def test_get_paps(self):
        paps = self.wordlist.get_paps(ref="cogid", loans=True)
        cogs = self.wordlist.get_etymdict(ref="cogid", loans=True)
        
        for key in cogs:
            if abs(key) in paps:
                assert True
            else:
                print(key)
                assert False
    
    def test_output(self):
        fn = text_type(self.tmp_path('test'))
        for fmt in 'taxa tre dst starling paps.nex paps.csv'.split():
            kw = {'ref': 'word'} if fmt == 'starling' else {}
            self.wordlist.output(fmt, filename=fn, **kw)
