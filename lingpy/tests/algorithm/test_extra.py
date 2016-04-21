from __future__ import unicode_literals
from lingpy.algorithm.extra import *
from mock import MagicMock, patch
from collections import namedtuple

class components():
    def __init__(self, nodes):
        self.vs = nodes
    def subgraphs(self):
        return [Igraph.Graph([v['name']]) for v in self.vs]

class Igraph(MagicMock):
    pass
    class Graph():
        def __init__(self, vs=[]):
            self.vs = [{'name':v} for v in vs]
        def add_edge(self, a, b):
            pass
        def add_vertex(self, vertex):
            self.vs += [{'name': vertex}]
        def community_infomap(self, *args, **kw):
            return components(self.vs)
class Cluster(MagicMock):
    pass
    def dbscan(self, *args, **kw):
        return (None, [-1 for i in range(len(args[0]))])
    class AffinityPropagation():
        def __init__(self, *args, **kw):
            pass
        def fit_predict(self, arg):
            return [i for i in range(len(arg))]


class Tests(object):
    def setUp(self):

        self.matrix = [[0.0, 0.5, 0.67, 0.8, 0.2],
                [0.5, 0.0, 0.4, 0.7, 0.6],
                [0.67, 0.4, 0.0, 0.8, 0.8],
                [0.8, 0.7, 0.8, 0.0, 0.3],
                [0.2, 0.6, 0.8, 0.3, 0.0]]
        self.taxa = ['German','Swedish','Icelandic','English','Dutch']
        
    @patch("lingpy.algorithm.extra.cluster", new=Cluster()) 
    @patch("lingpy.algorithm.extra.igraph", new=Igraph())
    def test_clustering(self):
        if not cluster:
            cluster1 = dbscan(0.25, self.matrix, self.taxa)
            cluster2 = dbscan(0.25, self.matrix, self.taxa, revert=True)
            cluster3 = affinity_propagation(0.5, self.matrix, self.taxa)
            cluster4 = affinity_propagation(0.5, self.matrix, self.taxa, revert=True)
            assert cluster2[0] != cluster2[4]
            assert cluster4[0] != cluster4[4]

        if not igraph:
            cluster5 = infomap_clustering(0.4, self.matrix, self.taxa)
            cluster6 = infomap_clustering(0.4, self.matrix, self.taxa, revert=True)
            assert cluster6[0] != cluster6[4]

    def test_dbscan(self):
        if cluster:
            cluster1 = dbscan(0.25, self.matrix, self.taxa)
            cluster2 = dbscan(0.25, self.matrix, self.taxa, revert=True)
            assert cluster2[0] == cluster2[4]

    def test_affinity_propagation(self):
        if cluster:
            cluster1 = affinity_propagation(0.5, self.matrix, self.taxa)
            cluster2 = affinity_propagation(0.5, self.matrix, self.taxa, revert=True)
            assert cluster2[0] == cluster2[4]
    
    def test_infomap_clustering(self):
        if igraph:
            cluster1 = infomap_clustering(0.4, self.matrix, self.taxa)
            cluster2 = infomap_clustering(0.4, self.matrix, self.taxa, revert=True)
            assert cluster2[0] == cluster2[4]

