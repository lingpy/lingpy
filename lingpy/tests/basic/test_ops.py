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

    def test_wl2qlc(self):
        from lingpy.basic.ops import wl2qlc

        stamp = 'test-stamp'
        out = self.tmp_path('test')
        wl2qlc(self.wordlist.header, self.wordlist._data, filename=str(out), stamp=stamp)
        out = self.tmp_path('test.qlc')
        with out.open(encoding='utf8') as fp:
            self.assertTrue(fp.read().endswith(stamp))

    def test_tsv2triple(self):
        from lingpy.basic.ops import tsv2triple, triple2tsv

        out = self.tmp_path('test')
        tsv2triple(self.wordlist, str(out))
        table = triple2tsv(str(out))
        self.assertIsInstance(table, list)

    def test_calculate_data(self):
        from lingpy.basic.ops import calculate_data

        for data in ['tree', 'dst', 'cluster']:
            calculate_data(self.wordlist, data)
            calculate_data(self.wordlist, data, mode='shared')
