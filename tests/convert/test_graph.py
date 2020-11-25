import lingpy.convert.graph as lcg

import networkx as nx

try:  # pragma: no cover
    import igraph as ig

    noigraph = False
except ImportError:  # pragma: no cover
    noigraph = True


def test_networkx2igraph():
    graph = nx.Graph()
    graph.add_nodes_from(['a', 'b', 'c'])
    graph.add_edges_from([('a', 'b'), ('c', 'd')])

    if not noigraph:  # pragma: no cover
        network = lcg.networkx2igraph(graph)
        assert len(network.vs) == 4
        assert len(network.es) == 2


def test_igraph2networkx():
    if not noigraph:  # pragma: no cover
        graph = ig.Graph()
        graph.add_vertices(['a', 'b', 'c', 'd'])
        graph.add_edges([('a', 'b'), ('c', 'd')])

        network = lcg.igraph2networkx(graph)
        assert len(list(network.nodes())) == 4
        assert len(list(network.edges())) == 2
