# created: Mo 21 Jan 2013 04:05:57  CET
# modified: Mo 21 Jan 2013 04:05:57  CET

__author__ = "Johann-Mattis List"
__date__ = "2013-01-21"

"""
Conversion routines for the GML format.
"""

from ..check.exceptions import ThirdPartyModuleError
try:
    import networkx as nx
except:
    ThirdPartyModuleError('networkx').warning()

from ..thirdparty import cogent as cg
from ..check.messages import FileWriteMessage

def gls2gml(
        gls,
        graph,
        tree,
        filename = 'gml'
        ):
    """
    Create GML-representation of a given gain-loss-scenario (GLS).

    Parameters
    ----------
    gls : list
        A list of tuples, indicating the origins of characters along a tree.
    graph : networkx.graph
        A graph that serves as a template for the plotting of the GLS.
    tree : cogent.tree.PhyloNode
        A tree object. 
    """
    
    # create a mapper for the ids and the string-names
    mapper = {}
    for node,data in graph.nodes(data=True):
        mapper[data['label']] = node

    # create a graph
    g = nx.Graph()

    # sort the gls according to the number of tips
    gls_srt = sorted(
            gls,
            key=lambda x:len(tree.getNodeMatchingName(x[0]).tips()),
            reverse=True
            )

    # set the basic event frame, depending on the state of the root
    if gls_srt[0][1] == 1 and gls_srt[0][0] == 'root':
        this_color = "#ffffff"
    else:
        this_color = "#000000"

    # let all nodes inherit these parameters
    for node,data in graph.nodes(data=True):
        data['graphics']['fill'] = this_color
        data['graphics']['type'] = 'ellipse'
        data['graphics']['w'] = 20.0
        data['graphics']['h'] = 20.0
        data['origin'] = 0
        g.add_node(node,**data)

    # assign the root as starting point
    data = graph.node[mapper['root']]
    data['graphics']['type'] = 'ellipse'
    data['graphics']['w'] = 50.0
    data['graphics']['h'] = 50.0
    g.add_node(mapper['root'],**data)

    # iterate over the nodes involved in change and assign the values to their
    # children
    for name,event in gls_srt:
        if event == 1:
            this_fill = '#ffffff'
        else:
            this_fill = '#000000'

        # get the names of the descendant nodes in the subtree 
        sub_tree_nodes = tree.getNodeMatchingName(name).getNodeNames()

        # iterate over all nodes to change
        for node in sub_tree_nodes:
            data = g.node[mapper[node]]
            data['graphics']['fill'] = this_fill
            g.add_node(mapper[node],**data)

        # change the size of the root of the subtree
        g.node[mapper[name]]['graphics']['h'] = 50.0
        g.node[mapper[name]]['graphics']['w'] = 50.0
        g.node[mapper[name]]['graphics']['fill'] = this_fill
        g.node[mapper[name]]['origin'] = 1

    # add the edges to the tree
    for edgeA,edgeB,data in graph.edges(data=True):
        # for computers with new networkx version
        try:
            del data['graphics']['Line']
        except:
            pass
        #if 'label' not in data:
        g.add_edge(edgeA,edgeB,**data)
    
    f = open(filename+'.gml','w')
    for line in nx.generate_gml(g):
        f.write(line+'\n')
    f.close()

    FileWriteMessage(filename,'gml').message('written')
    
    return g

def nwk2gml(
        treefile,
        filename = 'nwk'
        ):
    """
    Function converts a tree in newick format to a network in gml-format.
    """
    
    # create an empty graph
    graph = nx.DiGraph()
    
    # load the tree
    tree = cg.LoadTree(treefile)

    # get the node names of the tree
    nodes = tree.getNodeNames()

    # iterate over the nodes and add them and the edges to the graph
    for node in nodes:
        
        # add the node (just as a precaution)
        graph.add_node(node)

        # get the parent of the node
        parent = tree.getNodeMatchingName(node).Parent

        # add the edge if the parent is not None
        if parent:
            graph.add_edge(parent.Name,node)
    
    f = open(filename+'.gml','w')
    for line in nx.generate_gml(graph):
        f.write(line+'\n')
    f.close()

    FileWriteMessage(filename,'gml').message('written')

    return 
