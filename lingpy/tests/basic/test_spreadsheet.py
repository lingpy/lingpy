from unittest import TestCase

from lingpy.tests.util import test_data

# check for regex, otherwise no use to run this test
try:
    import regex
    with_regex = True
except ImportError:
    with_regex = False

class TestSpreadsheet(TestCase):
    def setUp(self):
        from lingpy.basic.spreadsheet import Spreadsheet

        self.spreadsheet = Spreadsheet(test_data('leach1969-67-161.csv'))

    def test_analyze(self):
        if with_regex:
            assert self.spreadsheet.analyze('characters', 'graphemes', 'words')
