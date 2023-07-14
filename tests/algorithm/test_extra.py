import pytest

from lingpy.algorithm.extra import *


@pytest.fixture
def Igraph(mocker):
    class Components:
        def __init__(self, nodes):
            self.vs = nodes

        def subgraphs(self):
            return [Igraph_.Graph([v['name']]) for v in self.vs]

    class Igraph_(mocker.MagicMock):
        pass

        class Graph:
            def __init__(self, vs=[]):
                self.vs = [{'name': str(v)} for v in vs]

            def add_edge(self, a, b):
                pass

            def add_vertex(self, vertex):
                self.vs += [{'name': str(vertex)}]

            def community_infomap(self, *args, **kw):
                return Components(self.vs)

    return Igraph_()


@pytest.fixture
def Cluster(mocker):
    class Cluster_(mocker.MagicMock):
        @staticmethod
        def dbscan(*args, **kw):
            return None, [-1 for _ in range(len(args[0]))]

        class AffinityPropagation:
            def __init__(self, *args, **kw):
                pass

            @staticmethod
            def fit_predict(arg):
                return [i for i in range(len(arg))]

    return Cluster_()


@pytest.fixture
def matrix():
    return [[0.0, 0.5, 0.67, 0.8, 0.2],
            [0.5, 0.0, 0.4, 0.7, 0.6],
            [0.67, 0.4, 0.0, 0.8, 0.8],
            [0.8, 0.7, 0.8, 0.0, 0.3],
            [0.2, 0.6, 0.8, 0.3, 0.0]]


@pytest.fixture
def taxa():
    return ['German', 'Swedish', 'Icelandic', 'English', 'Dutch']


def test_clustering(mocker, Cluster, Igraph, matrix, taxa):
    mocker.patch("lingpy.algorithm.extra.cluster", new=Cluster)
    mocker.patch("lingpy.algorithm.extra.igraph", new=Igraph)

    if not cluster:
        cluster1 = dbscan(0.25, matrix, taxa, revert=True)
        cluster2 = affinity_propagation(0.5, matrix, taxa, revert=True)
        assert cluster1[0] != cluster1[4]
        assert cluster2[0] != cluster2[4]

    if not igraph:
        cluster3 = infomap_clustering(0.4, matrix, taxa, revert=True)
        assert cluster3['0'] != cluster3['4']


def test_dbscan(matrix, taxa):
    if cluster:  # pragma: no cover
        cluster1 = dbscan(0.25, matrix, taxa, revert=True)
        assert cluster1[0] == cluster1[4]


def test_affinity_propagation(matrix, taxa):
    if cluster:  # pragma: no cover
        cluster1 = affinity_propagation(0.5, matrix, taxa, revert=True)
        assert cluster1[0] == cluster1[4]


def test_infomap_clustering(matrix, taxa):
    if igraph:  # pragma: no cover
        cluster1 = infomap_clustering(0.4, matrix, taxa, revert=True)
        assert cluster1['0'] == cluster1['4']
