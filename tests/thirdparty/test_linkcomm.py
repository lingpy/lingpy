from lingpy.thirdparty.linkcomm import link_clustering


def test_hlc():
    hlc = link_clustering.HLC(
        {"a": {"b", "c"},
         "b": {"c"},
         "c": {"a"}
         },
        {("a", "b"), ("b", "c"), ("a", "c")}
    )
    hlc.initialize_edges()
    assert hlc.cid2edges

    hlc.merge_comms(("a", "b"), ("b", "c"))
    hlc.single_linkage()
    tmp = hlc.single_linkage(threshold=0.5)[0]
    assert tmp['a', 'b'] == tmp['b', 'c']
