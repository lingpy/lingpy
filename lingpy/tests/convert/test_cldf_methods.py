from unittest import TestCase

from lingpy.basic.wordlist import Wordlist
from lingpy.tests.util import test_data


class Tests(TestCase):
    def test_from_cldf(self):
        wl = Wordlist.from_cldf(
            test_data('cldf/test-metadata.json'),
            col="Language_ID".lower(),
            row="Parameter_ID".lower())

        assert wl.width == 29
        assert wl.height == 1
        assert wl.entries[0] == 'alignment'
        assert wl.cols[0] == 'anuta'.lower()
        assert wl.cols[28] == 'wallisian'
