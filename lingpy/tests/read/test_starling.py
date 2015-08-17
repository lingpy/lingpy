# coding: utf8
from __future__ import unicode_literals

from lingpy.tests.util import WithTempDir, test_data


class Tests(WithTempDir):
    def test_star2qlc(self):
        from lingpy.read.starling import star2qlc

        star2qlc(test_data('rom.starling.tsv'), debug=True)
        res = star2qlc(test_data('rom.starling.tsv'))
        self.assertEquals(len(res), 208)
