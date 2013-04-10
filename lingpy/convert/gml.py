# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-09 08:39
# modified : 2013-04-09 08:39
"""
Conversion routines for the GML format.
"""

__author__="Johann-Mattis List"
__date__="2013-04-09"

from ..check.exceptions import ThirdPartyModuleError
import numpy as np
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
        filename = 'lingpy',
        ):
    """
    Function converts a tree in newick format to a network in gml-format.

    treefile : str
        Either a str defining the path to a file containing the tree in
        Newick-format, or the tree-string itself.
    filename : str (default='lingpy')
        The name of the output GML-file. If filename is set to c{None}, the
        function returns a :py:class:`~networkx.Graph`.

    Returns
    -------
    graph : networkx.Graph

    """
    
    # create an empty graph
    graph = nx.DiGraph()
    
    # load the tree
    try:
        tree = cg.LoadTree(treefile)
    except:
        tree = cg.LoadTree(treestring=treefile)

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
    
    if filename:
        f = open(filename+'.gml','w')
        for line in nx.generate_gml(graph):
            f.write(line+'\n')
        f.close()

        FileWriteMessage(filename,'gml').message('written')
        return
    else:
        return graph

def radial_layout(
        treestring,
        change = lambda x:x**1.75,
        degree = 100,
        filename = ''
        ):
    """
    Function calculates a simple radial tree layout.

    Parameters
    ----------
    treefile : str
        Either a str defining the path to a file containing the tree in
        Newick-format, or the tree-string itself.
    filename : str (default=None)
        The name of the output file (GML-format). If set to c{None}, no output
        will be written to file.
    change : function (default = lambda x:2 * x**2)
        The function used to modify the radius in the polar projection of the
        tree.

    Returns
    -------
    graph : networkx.Graph
        A graph representation of the tree with coordinates specified in the
        graphics-attribute of the nodes.

    Note
    ----
    This function creates a radial tree-layout from a given tree specified in
    Newick format.

    """
    # calculate the factor for projection from the degree
    pfactor = degree / 360

    # calculate the projection (should be centered)
    if degree < 180:
        pstart = (180 - degree) / 360 * np.pi
        pend = pstart + 2 * np.pi * pfactor
    else:
        pstart = 0
        pend = 2 * np.pi * pfactor

    # define private function for centering of nodes
    def get_center(nodes):
        
        # first sort all values since we need max and min of the theta values
        xvals = sorted([n[0] for n in nodes])
        
        # get minimum and maximum
        xA,xB = xvals[0],xvals[-1]

        # calculate the new coordinates, the radius is simply decreased by 1
        y = min([n[1] for n in nodes]) - 1
        
        # the theta-value is calculated by the following formula
        x = ( xA + abs(xA - xB) / 2 )

        return x,y

    # get the tree
    try:
        tree = cg.LoadTree(treestring)
    except:
        tree = cg.LoadTree(treestring=treestring)

    # get the leaves
    leaves = tree.getTipNames()

    # get the paths in order to find out the radius of the tree
    paths = {}
    
    for l in leaves:
        path = tree.getConnectingEdges('root',l)
        try:
            paths[len(path)] += [l]
        except:
            paths[len(path)] = [l]

    # get the max path
    maxL = max(paths)

    # get the initial coordinates
    coords = {}

    for node,x in zip(leaves,np.linspace(pstart,pend,len(leaves)+1)):
        coords[node] = (x,maxL)

    # assign leaves to queue
    queue = [l for l in leaves]

    # make the visited list
    visited = []

    # start the loop
    while queue:
        
        # get the node
        node = queue.pop(0)

        if node in visited:
            pass
        else:

            # get the parent and all children
            children = [child.Name for child in
                    tree.getNodeMatchingName(node).Parent.Children]

            # iterate over children
            goon = True
            for child in children:
                if child in coords:
                    pass
                else:
                    goon = False
                    break

            # goon, if this is possible
            if not goon:
                queue += [node]
            else:

                x,y = get_center(
                       [coords[child] for child in children]
                       )
                parent = tree.getNodeMatchingName(node).Parent.Name
                coords[parent] = (x,y)

                visited += [child for child in children]
                
                if parent != 'root':
                    queue += [parent]
    
    # convert tree to graph
    graph = nwk2gml(treestring,filename=None)
    
    # iterate over the graph and assign the data
    for n,d in graph.nodes(data=True):
        x,y = coords[n]

        # change coordinates
        xN = change(y) * np.cos(x)
        yN = change(y) * np.sin(x)

        # get angle for text-rotation in degrees
        angle = x * 180 / np.pi

        # check for specific parts where the angle has to be adapted
        if 270 > angle > 90:
            angle += 180

        # assign the data to the graph
        d['graphics'] = {'x':xN,'y':yN,'angle':angle}

        # don't forget the label
        d['label'] = n
    
    #for node,(x,y) in coords.items():
    #    
    #    xN = change(y) * np.cos(x)
    #    yN = change(y) * np.sin(x)

    #    coords[node] = (xN,yN)

    #
    ## append coords to graph
    #for node,data in graph.nodes(data=True):
    #    x,y = coords[node]
    #    data['graphics'] = {'x':x,'y':y}
    #    data['label'] = node

    #    graph.add_node(node,**data)
    
    if filename:
        f = open(filename+'.gml','w')
        for line in nx.generate_gml(graph):
            f.write(line+'\n')
        f.close()
        FileWriteMessage(filename,'gml').message('written')

    return graph

