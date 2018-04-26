import pytest

from lingpy.algorithm._tree import _TreeDist as tD


def test_grf():
    tree_a = '((a,b),(c,d));'
    tree_b = '((a:1,c:1):1,(b:1,d:1):1);'
    tree_c = '(((a,b),c),((d,e),f));'
    tree_d = '((a,b),(c,(d,(e,f))));'

    assert tD.grf(tree_a, tree_b) == 1.0
    assert tD.grf(tree_c, tree_d) == pytest.approx(0.33, 0.1)

    with pytest.raises(ValueError):
        tD.grf(tree_a, '((a,b),e);')


def test_get_bipartition():
    _, _ = tD.get_bipartition('(a,b,c),(c,d)')
    _, _ = tD.get_bipartition('((a,b,c),(c,d))')

    with pytest.raises(ValueError):
        tD.get_bipartition('a,b')

    with pytest.raises(ValueError):
        tD.get_bipartition('((),(a,b))')
