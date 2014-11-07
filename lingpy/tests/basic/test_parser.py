from unittest import TestCase

from lingpy.tests.util import test_data


class TestParser(TestCase):
    def test_Parser(self):
        from lingpy.basic.parser import QLCParser

        parser = QLCParser(test_data('KSL.qlc'))
        assert len(parser)
