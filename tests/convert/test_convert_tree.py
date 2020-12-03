from lingpy.basic.tree import Tree
from lingpy.convert import tree


def test__nwk_format():
    assert tree._nwk_format("test (taxon!?)") == "test_taxon"


def test_nwk2tree_matrix():
    newick = '(((a,b),(c,d)),e);'
    matrix, taxa = tree.nwk2tree_matrix(newick)
    assert taxa == Tree(newick).taxa
