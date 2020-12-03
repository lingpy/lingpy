import pytest

from lingpy.basic.tree import Tree, random_tree, _star_tree


def test_star_tree():
    assert _star_tree(list('abc')) == '(a,b,c);'


def test_random_tree():
    taxa = list('abcdefg')

    t1 = Tree(random_tree(taxa, branch_lengths=False))
    t2 = Tree(random_tree(taxa, branch_lengths=True))

    assert sorted(t1.taxa) == sorted(t2.taxa)
    assert str(t1) != str(t2)
    assert ':' in str(t2) and ':' not in str(t1)


@pytest.fixture
def tree():
    return Tree('((a,b),(c,d),e)')


def test_basic(tree):
    assert sorted(tree.taxa) == list('abcde')


def test_getDistanceToRoot(tree):
    assert tree.getDistanceToRoot('a') == 2


def test_get_LCA(tree):
    assert str(tree.get_LCA('a', 'b')) == '(a,b);'


def test_get_distance():
    tree_a = Tree('((a:1,b:1):1,(c:1,d:1):1)')
    tree_b = Tree('((a:1,c:1):1,(b:1,d:1):1)')

    assert tree_a.get_distance(tree_b, 'grf') == 1.0
    assert tree_a.get_distance(tree_b, 'rf') == 1.0
    assert tree_a.get_distance(tree_b, 'symmetric') == 2
    assert tree_a.get_distance(tree_b, 'grf') == 1.0


def test_get_distance_unknown():
    """test failure with unknown distance"""
    with pytest.raises(ValueError):
        Tree('(a,b)').get_distance(Tree('(a,b)'), 'xxx')


def test_init_from_file(test_data):
    tree = Tree(str(test_data / 'phybo.tre'))
    assert len(tree.taxa) == 40, "should have a taxa attribute and 40 tips"


def test_init_from_list():
    tree = Tree(['Simon', 'Mattis', 'Robert'])
    assert len(tree.taxa) == 3, "should have a taxa attribute and 3 tips"
