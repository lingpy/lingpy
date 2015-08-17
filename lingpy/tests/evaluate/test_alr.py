# coding: utf8
from __future__ import unicode_literals

from lingpy import Wordlist
from lingpy.tests.util import test_data, WithTempDir


class Tests(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.wl = Wordlist(test_data('KSL.qlc'))

    def test_med(self):
        from lingpy.evaluate.alr import med

        self.assertAlmostEquals(
            med(self.wl, gold='gloss', test='gloss', classes=True), 0.0)
