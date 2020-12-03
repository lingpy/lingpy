import pytest

from lingpy.algorithm._tree import _TreeDist as td


def test_grf():
    treeA = '((a,b),(c,d));'
    treeB = '((a:1,c:1):1,(b:1,d:1):1);'
    treeC = '(((a,b),c),((d,e),f));'
    treeD = '((a,b),(c,(d,(e,f))));'

    assert td.grf(treeA, treeB) == 1.0
    with pytest.raises(ValueError):
        td.grf(treeA, '((a,b),e);')
    a = td.grf(treeC, treeD)
    assert '{0:.2}'.format(a) == '0.33'
    with pytest.raises(ValueError):
        td.get_bipartition('a,b')
    td.get_bipartition('(a,b,c),(c,d)')
    td.get_bipartition('((a,b,c),(c,d))')
    with pytest.raises(ValueError):
        td.get_bipartition('((),(a,b))')
