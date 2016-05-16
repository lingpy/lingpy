"""
Test wordlist module.
"""
from six import text_type
from nose.tools import assert_raises

from lingpy import Wordlist
from lingpy.tests.util import test_data, WithTempDir


class TestWordlist(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.wordlist = Wordlist(test_data('KSL.qlc'))
        self.wordlist2 = Wordlist(test_data('good_file.tsv'))

    def test___len__(self):
        assert len(self.wordlist) == 1400

    def test_calculate(self):
        self.wordlist.calculate('dst')

        assert hasattr(self.wordlist, 'distances')
        assert sum([self.wordlist.distances[x][x] for x in
                    range(self.wordlist.width)]) == 0

        self.wordlist.calculate('tree')
        assert str(self.wordlist.tree).endswith(';')

        assert sorted(self.wordlist.tree.taxa) == sorted(self.wordlist.cols)

        self.wordlist.calculate('groups')

        assert hasattr(self.wordlist, 'groups')
        assert type(self.wordlist.groups) == dict

    def test_coverage(self):
        self.wordlist.coverage()
        self.wordlist.coverage(stats='ratio')
        self.wordlist.coverage(stats='mean')

    def test_get_list(self):
        gerL = self.wordlist.get_list(doculect='German', entry='ipa', flat=True)
        gerD = self.wordlist.get_dict(col='German', entry='ipa')
        gerT = self.wordlist.get_list(doculect='German', entry="ipa")

        assert sorted(gerL) == sorted([v[0] for v in gerD.values()])
        assert sorted(gerT) == sorted(gerL)

        hand1 = self.wordlist.get_list(concept="hand", entry="ipa", flat=True)
        hand2 = self.wordlist.get_dict(row="hand", entry="ipa")
        hand3 = self.wordlist.get_list(concept="hand", flat=True)
        assert sorted(hand1) == sorted([v[0] for v in hand2.values()])

        # test for synonym lines, which are flattened
        assert self.wordlist2.get_list(concept='hand', entry="language",
                flat=True).count('l6') == 2
        nonflat = self.wordlist2.get_list(concept="hand", entry="language")
        assert nonflat[0][-1] == nonflat[1][-1]
        assert len(self.wordlist2.get_list(col="l1", entry="concept")) == 3
        assert len(self.wordlist2.get_list(col="l1", flat=True, entry="concept")) == 2

        assert_raises(ValueError, self.wordlist2.get_list, col="l1",
                row="hand")
        assert_raises(ValueError, self.wordlist2.get_list)
        assert_raises(ValueError, self.wordlist.get_list, **{"row" : "Hand"})

    def test_get_dict(self):
        gerD = self.wordlist.get_dict(col='German')

        assert sorted(gerD.keys()) == sorted(self.wordlist.rows)
        assert_raises(ValueError, self.wordlist.get_dict, **{"row" : "Hand"})

    def test_renumber(self):
        self.wordlist.renumber('cogid', 'dummy')

        ger1 = self.wordlist.get_list(col='German', entry='cogid', flat=True)
        ger2 = self.wordlist.get_list(col='German', entry='dummy', flat=True)

        assert len(set(ger1)) == len(set(ger2))
        assert sum([1 for x in ger2 if type(x) == int]) == len(ger2)

    def test_get_entries(self):
        ger = self.wordlist.get_entries('cogid')

        assert len(ger) == self.wordlist.height
        assert len(ger[0]) == self.wordlist.width

    def test_get_etymdict(self):
        etd1 = self.wordlist.get_etymdict(ref='cogid', entry='ipa', modify_ref=False)
        etd2 = self.wordlist.get_etymdict(ref='cogid', entry='ipa',
                modify_ref=abs)

        assert len(etd1) > len(etd2) and len(set([abs(x) for x in etd1])) == \
                                         len(etd2)
        assert len([x for x in etd2 if x < 0]) == 0

        # make "fuzzy" cognate sets
        self.wordlist.add_entries('fuzzyid', 'cogid', lambda x: [x])

        etd3 = self.wordlist.get_etymdict(
            ref='fuzzyid', entry='ipa', modify_ref=False)
        etd4 = self.wordlist.get_etymdict(
            ref='fuzzyid', entry='ipa', modify_ref=abs)
        for key in etd1:
            assert etd1[key] == etd3[key]
        for key in etd2:
            self.assertEquals(etd2[key], etd4[key])

    def test_get_paps(self):
        paps = self.wordlist.get_paps(ref="cogid", modify_ref=abs)
        cogs = self.wordlist.get_etymdict(ref="cogid", modify_ref=abs)

        for key in cogs:
            if abs(key) in paps:
                assert True
            else:
                print(key)
                assert False

    def test_output(self):
        fn = text_type(self.tmp_path('test'))
        for fmt in 'tsv taxa tre dst starling paps.nex paps.csv separated multistate.nex groups'.split():
            kw = {'ref': 'word'} if fmt == 'starling' else {}
            self.wordlist.output(fmt, filename=fn, **kw)
            if fmt == 'starling':
                self.wordlist.output(fmt, filename=fn, cognates='cogid', **kw)
            if fmt == 'tsv':
                kw['subset'] = True
                self.wordlist.output(fmt, filename=fn, cols=[], rows={}, **kw)
                self.wordlist.output(fmt, filename=fn,
                        cols=sorted(self.wordlist.header)[:2], rows=dict(ID=" > 10"),
                            **kw)
    def test_export(self):
        fn = text_type(self.tmp_path('test'))
        for fmt in 'txt tex html'.split():
            self.wordlist.export(fmt, filename=fn)

    def test_get_wordlist(self):
        from lingpy.basic.wordlist import get_wordlist
        wl1 = get_wordlist(test_data('mycsvwordlist.csv'))
        wl2 = get_wordlist(test_data('mycsvwordlistwithoutids.csv'))
        assert wl1.height == wl2.height
        for k in wl1:
            assert wl1[k, 'concept'] == wl2[k, 'concept']
