from unittest import TestCase

from lingpy.basic.tree import Tree, random_tree


def test_random_tree():
    taxa = list('abcdefg')

    t1 = Tree(random_tree(taxa, branch_lengths=False))
    t2 = Tree(random_tree(taxa, branch_lengths=True))

    assert sorted(t1.taxa) == sorted(t2.taxa)
    assert str(t1) != str(t2)
    assert ':' in str(t2) and ':' not in str(t1)


class TestTree(TestCase):
    def setUp(self):
        self.tree = Tree('((a,b),(c,d),e)')
        assert sorted(self.tree.taxa) == list('abcde')

    def test_getDistanceToRoot(self):
        assert self.tree.getDistanceToRoot('a') == 2

    def test_get_LCA(self):
        assert str(self.tree.get_LCA('a', 'b')) == '(a,b);'

    def test_get_distance(self):
        treeA = Tree('((a:1,b:1):1,(c:1,d:1):1)')
        treeB = Tree('((a:1,c:1):1,(b:1,d:1):1)')

        assert treeA.get_distance(treeB, 'grf') == 1.0
        assert treeA.get_distance(treeB, 'rf') == 1.0
        assert treeA.get_distance(treeB, 'symmetric') == 2
        assert treeA.get_distance(treeB, 'grf') == 1.0
