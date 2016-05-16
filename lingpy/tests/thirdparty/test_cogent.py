# *-* coding: utf-8 *-*
"""
Test thirdparty modules.
"""
from __future__ import unicode_literals, print_function, division
from unittest import TestCase
from nose.tools import assert_raises
from collections import defaultdict

from six import PY3

from lingpy.thirdparty.cogent import LoadTree
from lingpy.thirdparty.cogent.newick import TreeParseError
from lingpy.tests.util import test_data


def test_LoadTree():
    
    # test to load a given tree-file
    tree = LoadTree(test_data('phylogeny.tre'))
    
    taxa = sorted(["Beijing", "Changsha", "Chengdu", "Fuzhou",
            "Guangzhou", "Guiyang", "Haerbin", "Haikou", "Hangzhou", "Hefei",
            "Huhehaote", "Jian\u2019ou", "Jinan", "Kunming", "Lanzhou",
            "Meixian", "Nanchang", "Nanjing", "Nanning", "Pingyao", "Qingdao",
            "Shanghai", "Shantou", "Shexian", "Suzhou", "Taibei", "Taiyuan",
            "Taoyuan", "Tianjin", "Tunxi", "Wenzhou", "Wuhan", "Wulumuqi",
            "Xi\u2019an", "Xiamen", "Xianggang", "Xiangtan", "Xining",
            "Yinchuan", "Zhengzhou"])
    
    for a, b in zip(sorted(tree.taxa), taxa):
        assert a == b
    
    tree = LoadTree("((((((((Taiyuan,Pingyao,Huhehaote),((((Xi’an,Xining,Zhengzhou),(Lanzhou,Yinchuan,Wulumuqi)),(((Tianjin,Jinan),Qingdao),Beijing,Haerbin)),(((Guiyang,Kunming),Chengdu,Wuhan),(Nanjing,Hefei)))),(Xiangtan,Changsha)),Nanchang),(Shexian,Tunxi)),((Shanghai,Suzhou,Hangzhou),Wenzhou)),(((Xianggang,Guangzhou),Nanning),(Meixian,Taoyuan))),((((Xiamen,Taibei),Shantou,Haikou),Fuzhou),Jian’ou));")

    for a, b in zip(sorted(tree.taxa), taxa):
        assert a == b


class TreeTests(TestCase):
    def test_Tree(self):
        from lingpy.thirdparty.cogent import TreeNode, TreeError

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
        if not PY3:
            tree.makeTreeArray()
        tree.prune()
        tree.getDistances()
        list(tree.traverse_recursive())
        self.assertEquals(tree.lowestCommonAncestor(['b', 'c']).Name, 'a')
        self.assertEquals(tree.separation(tree[0]), 1)
        for before in [True, False]:
            for after in [True, False]:
                list(tree.traverse(before, after))
        assert tree.childGroups()
        tree[0].childGroups()
        tree.compareBySubsets(tree.sorted())
        for p1 in [True, False]:
            for p2 in [True, False]:
                tree.asciiArt(p1, p2, defaultdict(lambda: 'a'))
        self.assertRaises(TreeError, tree.getSubTree, [])
        tree.get_LCA(tree.getNodeNames()[0])
        tree.getNewickRecursive()
    
    def test_more_trees(self):

        tree = LoadTree('((a:1,b:1):2,(c:3,d:4):5);')
        tree2 = LoadTree('((a,b),(c,d));')
        assert tree.sameShape(tree)
        tree.tipToTipDistances()
        tree.getMaxTipTipDistance()
        tree.maxTipTipDistance()
        tree.getSubTree(['a','b','c'])
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


        ntree = tree.copyRecursive()
        ntree2 = tree2.copy()
        assert ntree.sameShape(ntree2)
        assert 'a' in tree.getNodeNames()

        assert_raises(TreeParseError, LoadTree, '(a,b;');
        assert_raises(TreeParseError, LoadTree, treestring="(a,b'");
        assert_raises(TreeParseError, LoadTree, treestring="(a,'b'");
        assert_raises(TreeParseError, LoadTree, treestring="(a,'b'\n");
        assert_raises(TreeParseError, LoadTree, treestring="(a,'b-c'");
        assert_raises(TreeParseError, LoadTree, "(a,b),('c_d',e);")
        new_tree = LoadTree("((a,b),('',e));")
        new_tree.nameUnnamedNodes()


class PhyloNodeTests(TestCase):
    def test_PhyloNode(self):
        from lingpy.thirdparty.cogent import PhyloNode

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
        if not PY3:
            node.compareByPartitions(PhyloNode(Name='x'))
        node.getDistances()
        node.unrootedDeepcopy()
        node.bifurcating()
        node.balanced()
        node.sameTopology(node)
        node.getEdgeNames('a', 'b', False, False)
        node.getEdgeNames('d', 'b', False, False, outgroup_name='c')
        node.getEdgeNames('d', 'b', True, True, outgroup_name='c')
        node._getDistances(endpoints=["d", "b"])

        
