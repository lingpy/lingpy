from unittest import TestCase

from lingpy.tests.util import test_data


class TestTree(TestCase):
    def test_init(self):
        from lingpy.basic.tree import Tree

        Tree([])
