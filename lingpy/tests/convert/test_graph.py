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


def test_networkx2igraph():
    if not noigraph and not nonetworkx:
        G = nx.Graph()
        G.add_nodes_from(['a', 'b', 'c'])
        G.add_edges_from([('a', 'b'), ('c', 'd')])

        N = lcg.networkx2igraph(G)
        assert len(N.vs) == 4
        assert len(N.es) == 2


def test_igraph2networkx():
    if not noigraph and not nonetworkx:
        G = ig.Graph()
        G.add_vertices(['a', 'b', 'c', 'd'])
        G.add_edges([('a', 'b'), ('c', 'd')])

        N = lcg.igraph2networkx(G)
        assert len(list(N.nodes())) == 4
        assert len(list(N.edges())) == 2
