"""
Test thirdparty modules.
"""
from collections import defaultdict

import pytest

from lingpy.thirdparty.cogent import LoadTree, TreeNode, TreeError, PhyloNode
from lingpy.thirdparty.cogent.newick import TreeParseError


def test_load_tree(test_data):
    # test to load a given tree-file
    tree = LoadTree(str(test_data / 'phylogeny.tre'))
    assert tree.compareByPartitions(tree) == 0
    taxa = sorted(["Beijing", "Changsha", "Chengdu", "Fuzhou", "Guangzhou",
                   "Guiyang", "Haerbin", "Haikou", "Hangzhou", "Hefei",
                   "Huhehaote", "Jian\u2019ou", "Jinan", "Kunming",
                   "Lanzhou",
                   "Meixian", "Nanchang", "Nanjing", "Nanning", "Pingyao",
                   "Qingdao", "Shanghai", "Shantou", "Shexian", "Suzhou",
                   "Taibei", "Taiyuan", "Taoyuan", "Tianjin", "Tunxi",
                   "Wenzhou", "Wuhan", "Wulumuqi", "Xi\u2019an", "Xiamen",
                   "Xianggang", "Xiangtan", "Xining", "Yinchuan",
                   "Zhengzhou"])
    assert len(tree.getConnectingEdges('Beijing', 'Xiangtan')) == 8
    for a, b in zip(sorted(tree.taxa), taxa):
        assert a == b
    tree.removeNode('Beijing')
    tree.prune()

    tree = LoadTree(tip_names=['Beijing', 'Changsha'])
    assert len(tree.taxa) == 2

    tree = LoadTree("((((((((Taiyuan,Pingyao,Huhehaote),"
                    "((((Xi’an,Xining,Zhengzhou),(Lanzhou,Yinchuan,"
                    "Wulumuqi)),"
                    "(((Tianjin,Jinan),Qingdao),Beijing,Haerbin)),"
                    "(((Guiyang,Kunming),Chengdu,Wuhan),(Nanjing,Hefei)))),"
                    "(Xiangtan,Changsha)),Nanchang),(Shexian,Tunxi)),"
                    "((Shanghai,Suzhou,Hangzhou),Wenzhou)),"
                    "(((Xianggang,Guangzhou),Nanning),(Meixian,Taoyuan))),"
                    "((((Xiamen,Taibei),Shantou,Haikou),Fuzhou),Jian’ou));")

    for a, b in zip(sorted(tree.taxa), taxa):
        assert a == b


def test_Tree():
    tree = TreeNode(
        Name='a',
        Children=[
            TreeNode(Name='b', Children=[TreeNode(Name='c')]),
            TreeNode(Name='b2'),
        ])
    tree.getEdgeNames('b', 'c', False, False)
    assert tree[0]
    tree.compareByNames(tree[0])
    '%s' % tree.copy()
    assert tree.params == {}
    tree.prune()
    tree.getDistances()
    list(tree.traverse_recursive())
    assert tree.lowestCommonAncestor(['b', 'c']).Name == 'a'
    assert tree.separation(tree[0]) == 1
    for before in [True, False]:
        for after in [True, False]:
            list(tree.traverse(before, after))
    assert tree.childGroups()
    tree[0].childGroups()
    tree.compareBySubsets(tree.sorted())
    for p1 in [True, False]:
        for p2 in [True, False]:
            tree.asciiArt(p1, p2, defaultdict(lambda: 'a'))
    with pytest.raises(TreeError):
        tree.getSubTree([])
    tree.get_LCA(tree.getNodeNames()[0])
    tree.getNewickRecursive()
    tree[0] = tree[0]


def test_more_trees():
    tree = LoadTree('((a:1,b:1):2,(c:3,d:4):5);')
    tree2 = LoadTree('((a,b),(c,d));')

    assert tree.sameShape(tree)

    tree.tipToTipDistances()
    tree.getMaxTipTipDistance()
    tree.maxTipTipDistance()
    tree.getSubTree(['a', 'b', 'c'])

    assert tree.compareName(tree2) == 0
    assert tree.compareName(tree) == 0
    assert tree.compareName('(a,b),(c,d));') != 0

    tree.descendantArray()
    tree.nameUnnamedNodes()
    tree.makeTreeArray()
    tree.getNewickRecursive(with_distances=True, semicolon=False,
                            escape_name=False)

    assert 'a' in tree.getNodesDict()
    assert 'b' in tree.taxa
    assert tree.getDistances()['a', 'b'] == 2
    assert tree.getMaxTipTipDistance()[0] == 12
    assert tree.maxTipTipDistance()[0] == 12
    tree.scaleBranchLengths(ultrametric=True)
    assert tree.getDistances()['a', 'b'] == pytest.approx(22.0)

    ntree = tree.copyRecursive()
    ntree2 = tree2.copy()
    assert ntree.sameShape(ntree2)
    assert 'a' in tree.getNodeNames()

    with pytest.raises(TreeParseError):
        LoadTree('(a,b;')
    with pytest.raises(TreeParseError):
        LoadTree(treestring="(a,b'")
    with pytest.raises(TreeParseError):
        LoadTree(treestring="(a,'b'")
    with pytest.raises(TreeParseError):
        LoadTree(treestring="(a,'b'\n")
    with pytest.raises(TreeParseError):
        LoadTree(treestring="(a,'b-c'")
    with pytest.raises(TreeParseError):
        LoadTree("(a,b),('c_d',e);")
    new_tree = LoadTree("((a,b),('',e));")
    new_tree.nameUnnamedNodes()


def test_PhyloNode():
    node = PhyloNode(
        Length=1,
        Name='a',
        Children=[
            PhyloNode(Length=7, Name='b'),
            PhyloNode(Length=3, Name='c')])
    node.scaleBranchLengths()
    node.append(PhyloNode(Length=7, Name='d'))
    node.rootAtMidpoint()
    node.unrooted()
    node.getDistances()
    node.unrootedDeepcopy()
    node.bifurcating()
    node.balanced()
    node.sameTopology(node)
    node.getEdgeNames('a', 'b', False, False)
    node.getEdgeNames('d', 'b', False, False, outgroup_name='c')
    node.getEdgeNames('d', 'b', True, True, outgroup_name='c')
    node._getDistances(endpoints=["d", "b"])
