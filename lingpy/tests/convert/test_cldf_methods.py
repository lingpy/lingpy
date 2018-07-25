from unittest import TestCase

from lingpy.basic.wordlist import Wordlist
from lingpy.tests.util import test_data
from lingpy.tests.util_testing import WithTempDir


class Tests(TestCase):
    def test_load_from_cldf_metadata(self):
        wl = Wordlist.from_cldf(
            test_data('cldf/test-metadata.json'),
            col="Language_ID".lower(),
            row="Parameter_ID".lower())

        assert wl.width == 29
        assert wl.height == 1
        assert wl.entries[0] == 'alignment'
        assert wl.cols[0] == 'anuta'.lower()
        assert wl.cols[28] == 'wallisian'

class CLDFWordlistWriteTest(WithTempDir):
    def test_load_cldf_and_write(self):
        wl = Wordlist.from_cldf(
            test_data('cldf/test-metadata.json'),
            col="Language_ID".lower(),
            row="Parameter_ID".lower())
        wl.output('tsv', filename=str(self.tmp_path('lingpycldf')))
