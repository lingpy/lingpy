from six import text_type
from lingpy.tests.util import WithTempDir 
from lingpy.tests.util import test_data
from lingpy.convert import tree
from lingpy.basic.tree import Tree

class TestTree(WithTempDir):
    
    def setUp(self):
        WithTempDir.setUp(self)
        self.newick = '(((a,b),(c,d)),e);'
        self.tree = Tree(self.newick)

    def test__nwk_format(self):

        assert tree._nwk_format("test (taxon!?)") == "test_taxon"

    def test_nwk2tree_matrix(self):
        
        matrix, taxa = tree.nwk2tree_matrix(self.newick)
        assert taxa == self.tree.taxa

