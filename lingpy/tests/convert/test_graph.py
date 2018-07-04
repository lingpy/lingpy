from unittest import TestCase

import lingpy.convert.graph as lcg

try:
    import networkx as nx

    nonetworkx = False
except ImportError:
    nonetworkx = True

try:
    import igraph as ig

    noigraph = False
except ImportError:
    noigraph = True


class Tests(TestCase):
    def test_networkx2igraph(self):
        if not noigraph and not nonetworkx:
            graph = nx.Graph()
            graph.add_nodes_from(['a', 'b', 'c'])
            graph.add_edges_from([('a', 'b'), ('c', 'd')])

            network = lcg.networkx2igraph(graph)
            assert len(network.vs) == 4
            assert len(network.es) == 2

    def test_igraph2networkx(self):
        if not noigraph and not nonetworkx:
            graph = ig.Graph()
            graph.add_vertices(['a', 'b', 'c', 'd'])
            graph.add_edges([('a', 'b'), ('c', 'd')])

            network = lcg.igraph2networkx(graph)
            assert len(list(network.nodes())) == 4
            assert len(list(network.edges())) == 2
