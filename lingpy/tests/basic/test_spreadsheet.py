from unittest import TestCase

from lingpy.tests.util import test_data


class TestSpreadsheet(TestCase):
    def setUp(self):
        from lingpy.basic.spreadsheet import Spreadsheet

        self.spreadsheet = Spreadsheet(test_data('leach1969-67-161.csv'))

    def test_analyze(self):
        assert self.spreadsheet.analyze('characters', 'graphemes', 'words')
