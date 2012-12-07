"""
Various functions for network conversion.
"""
__author__="Johann-Mattis List"
__date__="2012-12-06"

import networkx as nx
from ..thirdparty import cogent as cg

def tre2gml(infile):

    # create a digraph
    G = nx.DiGraph()

    # load the tree
    tree = cg.LoadTree(infile)
    
    # get the nodes
    nodes = tree.getNodeNames()

    # iterate over the nodes and add them to the graph
    for node in nodes:
        
        # add the node
        G.add_node(node)
        
        # get the parent
        parent = tree.getNodeMatchingName(node).Parent

        # add the edge if parent is not none
        if parent:
            G.add_edge(parent.Name,node)

    return G
