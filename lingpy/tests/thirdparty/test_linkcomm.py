from lingpy.thirdparty.linkcomm import link_clustering
from nose.tools import assert_raises
from lingpy.tests.util import test_data, WithTempDir, get_log

def test_swap():

    assert link_clustering.swap('b', 'a') == ("a", "b")
    assert link_clustering.swap(1, 2) == (1, 2)

def test_dc():

   assert int(link_clustering.Dc(10,5)) == 5
   assert int(link_clustering.Dc(10, 1)) == 0

class TestHLC(WithTempDir):

    def setUp(self):
        WithTempDir.setUp(self)
        self.hlc = link_clustering.HLC(
                { "a" : set(["b", "c"]),
                    "b" : set(["c"]),
                    "c" : set(["a"])
                    },
                set([("a", "b"), ("b", "c"), ("a", "c")])
                )
    def test_hlc(self):

        self.hlc.initialize_edges()
        assert self.hlc.cid2edges

        self.hlc.merge_comms(("a", "b"), ("b", "c"))
        self.hlc.single_linkage()
        tmp = self.hlc.single_linkage(threshold=0.5)[0]
        assert tmp['a','b'] == tmp['b','c']
        


        
