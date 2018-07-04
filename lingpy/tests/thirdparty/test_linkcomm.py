from lingpy.tests.util_testing import WithTempDir
from lingpy.thirdparty.linkcomm import link_clustering


class Tests(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.hlc = link_clustering.HLC(
            {"a": {"b", "c"},
             "b": {"c"},
             "c": {"a"}
             },
            {("a", "b"), ("b", "c"), ("a", "c")}
        )

    def test_hlc(self):
        self.hlc.initialize_edges()
        assert self.hlc.cid2edges

        self.hlc.merge_comms(("a", "b"), ("b", "c"))
        self.hlc.single_linkage()
        tmp = self.hlc.single_linkage(threshold=0.5)[0]
        assert tmp['a', 'b'] == tmp['b', 'c']
