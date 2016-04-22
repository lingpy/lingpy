"""
Test wordlist module.
"""
from __future__ import unicode_literals, print_function, division, absolute_import

from lingpy import Wordlist
from lingpy.tests.util import test_data, WithTempDir


class TestOps(WithTempDir):

    def setUp(self):
        WithTempDir.setUp(self)
        self.wordlist = Wordlist(test_data('KSL.qlc'))

    def test_wl2dict(self):
        from lingpy.basic.ops import wl2dict

        res = wl2dict(self.wordlist, dict(s1=['concept','{0}'], s2=['cogid',
            '{0}']), [('ipa', '{0}')])

    def test_wl2dst(self):
        from lingpy.basic.ops import wl2dst

        res = wl2dst(self.wordlist, mode='jaccard')
        self.assertIsInstance(res, list)
        res = wl2dst(self.wordlist, mode='jaccard', refB='glossid')
        self.assertIsInstance(res, list)

    def test_wl2qlc(self):
        from lingpy.basic.ops import wl2qlc

        stamp = 'test-stamp'
        out = self.tmp_path('test')

        wl2qlc(self.wordlist.header, self.wordlist._data, filename=out.as_posix(), stamp=stamp)
        out = self.tmp_path('test.qlc')
        with out.open(encoding='utf8') as fp:
            self.assertTrue(fp.read().endswith(stamp))

    def test_tsv2triple(self):
        from lingpy.basic.ops import tsv2triple, triple2tsv

        out = self.tmp_path('test')
        triples = tsv2triple(self.wordlist, out)
        self.assertIsInstance(triple2tsv(out), list)
        self.assertIsInstance(triple2tsv(triples, output='dict'), dict)

    def test_calculate_data(self):
        from lingpy.basic.ops import calculate_data

        for data in ['tree', 'dst', 'cluster']:
            calculate_data(self.wordlist, data)
            calculate_data(self.wordlist, data, mode='shared')

    def test_wl2multistate(self):
        from lingpy.basic.ops import wl2multistate

        res = wl2multistate(self.wordlist, 'cogid')
        # the result must be a matrix.
        self.assertIsInstance(res, list)
        self.assertEquals(len(set(len(row) for row in res)), 1)

    def test_coverage(self):
        from lingpy.basic.ops import coverage

        res = coverage(self.wordlist)
        self.assertEquals(res['Turkish'], 200)
