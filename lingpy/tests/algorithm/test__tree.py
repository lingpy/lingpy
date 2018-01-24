from nose.tools import assert_raises
from lingpy.algorithm._tree import _TreeDist as td
class Tests():

    def setUp(self):
        self.treeA = '((a,b),(c,d));'
        self.treeB = '((a:1,c:1):1,(b:1,d:1):1);'
        self.treeC = '(((a,b),c),((d,e),f));'
        self.treeD = '((a,b),(c,(d,(e,f))));'

    def test_grf(self):
        assert td.grf(self.treeA, self.treeB) == 1.0
        assert_raises(ValueError, td.grf, self.treeA, '((a,b),e);')
        a = td.grf(self.treeC, self.treeD)
        assert '{0:.2}'.format(a) == '0.33'
        assert_raises(ValueError, td.get_bipartition, 'a,b')
        td.get_bipartition('(a,b,c),(c,d)')
        td.get_bipartition('((a,b,c),(c,d))')
        assert_raises(ValueError, td.get_bipartition, '((),(a,b))')

