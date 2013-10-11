# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-10-07 14:05
# modified : 2013-10-07 14:05
"""
Module provides general clustering functions for LingPy.
"""

__author__="Johann-Mattis List"
__date__="2013-10-07"

try:
    from .cython import cluster as cluster
    from .cython import misc as misc
except ImportError:
    from .cython import _cluster as cluster
    from .cython import _misc as misc

from ..thirdparty import linkcomm as lc
from ..thirdparty import cogent as cg

from ..settings import rcParams

# thirdparty modules
import numpy as np

try:
    import networkx as nx
except ImportError:
    print(rcParams['W_missing_module'].format('networkx'))

def flat_upgma(threshold,matrix,taxa=[],revert=False):
    """
    Carry out a flat cluster analysis based on the UPGMA algorithm \
    (:evobib:`Sokal1958`).
    
    Parameters
    ----------
 
    threshold : float
        The threshold which terminates the algorithm.   
    
    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    taxa : list (default = [])
        A list containing the names of the taxa. If the list is left empty, the
        indices of the taxa will be returned instead of their names.
    
    Returns
    -------
    
    clusters : dict
        A dictionary with cluster-IDs as keys and a list of the taxa
        corresponding to the respective ID as values.

    Examples
    --------
    The function is automatically imported along with LingPy.

    >>> from lingpy import *
    
    Create a list of arbitrary taxa.

    >>> taxa = ['German','Swedish','Icelandic','English','Dutch']
    
    Create an arbitrary distance matrix.

    >>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])
    >>> matrix
    array([[ 0.  ,  0.5 ,  0.67,  0.8 ,  0.2 ],
           [ 0.5 ,  0.  ,  0.4 ,  0.7 ,  0.6 ],
           [ 0.67,  0.4 ,  0.  ,  0.8 ,  0.8 ],
           [ 0.8 ,  0.7 ,  0.8 ,  0.  ,  0.3 ],
           [ 0.2 ,  0.6 ,  0.8 ,  0.3 ,  0.  ]]

    Carry out the flat cluster analysis.

    >>> flat_upgma(0.5,matrix,taxa)
    {0: ['German', 'Dutch', 'English'], 1: ['Swedish', 'Icelandic']}

    See also
    --------
    lingpy.algorithm.clusters.upgma
    lingpy.algorithm.clusters.neighbor

    """
    
    return cluster.flat_upgma(threshold,matrix,taxa,revert)

def flat_cluster(
        method,
        threshold,
        matrix,
        taxa=[],
        revert=False
        ):
    """
    Carry out a flat cluster analysis based on the UPGMA algorithm.
    
    Parameters
    ----------
    method : str { 'upgma', 'single', 'complete' }
        Select between 'ugpma', 'single', and 'complete'.
 
    threshold : float
        The threshold which terminates the algorithm.   
    
    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    taxa : list (default = [])
        A list containing the names of the taxa. If the list is left empty, the
        indices of the taxa will be returned instead of their names.
    
    Returns
    -------
    
    clusters : dict
        A dictionary with cluster-IDs as keys and a list of the taxa
        corresponding to the respective ID as values.

    Examples
    --------
    The function is automatically imported along with LingPy.

    >>> from lingpy import *
    
    Create a list of arbitrary taxa.

    >>> taxa = ['German','Swedish','Icelandic','English','Dutch']
    
    Create an arbitrary distance matrix.

    >>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])
    >>> matrix
    array([[ 0.  ,  0.5 ,  0.67,  0.8 ,  0.2 ],
           [ 0.5 ,  0.  ,  0.4 ,  0.7 ,  0.6 ],
           [ 0.67,  0.4 ,  0.  ,  0.8 ,  0.8 ],
           [ 0.8 ,  0.7 ,  0.8 ,  0.  ,  0.3 ],
           [ 0.2 ,  0.6 ,  0.8 ,  0.3 ,  0.  ]])

    Carry out the flat cluster analysis.

    >>> flat_upgma(0.5,matrix,taxa)
    {0: ['German', 'Dutch', 'English'], 1: ['Swedish', 'Icelandic']}

    See also
    --------
    lingpy.algorithm.clusters.upgma
    lingpy.algorithm.clusters.neighbor

    """
    return cluster.flat_cluster(method,threshold,matrix,taxa,revert)

def upgma(
        matrix,
        taxa,
        distances = True
        ):
    """
    Carry out a cluster analysis based on the UPGMA algorithm \
    (:evobib:`Sokal1958`).

    Parameters
    ----------

    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    taxa : list
        An list containing the names of all taxa corresponding to the distances
        in the matrix.

    distances : bool
        If set to ``False``, only the topology of the tree will be returned.

    Returns
    -------

    newick : str
        A string in newick-format which can be further used in biological
        software packages to view and plot the tree.

    Examples
    --------
    Function is automatically imported when importing lingpy.

    >>> from lingpy import *
    
    Create an arbitrary list of taxa.

    >>> taxa = ['German','Swedish','Icelandic','English','Dutch']
    
    Create an arbitrary matrix.

    >>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])

    Carry out the cluster analysis.

    >>> upgma(matrix,taxa,distances=False)
    '((Swedish,Icelandic),(English,(German,Dutch)));'

    See also
    --------
    lingpy.algorithm.cluster.neighbor
    lingpy.algorithm.cluster.flat_upgma
   
    """

    return cluster.upgma(matrix,taxa,distances)

def neighbor(matrix,taxa,distances=True):
    """
    Function clusters data according to the Neighbor-Joining algorithm \
    (:evobib:`Saitou1987`).
    
    Parameters
    ----------

    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    taxa : list
        An list containing the names of all taxa corresponding to the distances
        in the matrix.

    distances : bool
        If set to ``False``, only the topology of the tree will be returned.

    Returns
    -------

    newick : str
        A string in newick-format which can be further used in biological
        software packages to view and plot the tree.

    Examples
    --------
    Function is automatically imported when importing lingpy.

    >>> from lingpy import *
    
    Create an arbitrary list of taxa.

    >>> taxa = ['Norwegian','Swedish','Icelandic','Dutch','English']
    
    Create an arbitrary matrix.

    >>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])

    Carry out the cluster analysis.

    >>> neighbor(matrix,taxa)
    '(((Norwegian,(Swedish,Icelandic)),English),Dutch);'

    See also
    --------
    lingpy.algorithm.cluster.upgma
    lingpy.algorithm.cluster.flat_upgma
    """
    
    return cluster.upgma(matrix,taxa,distances)

def fuzzy(threshold,matrix,taxa,method='upgma',revert=False):
    """
    Create fuzzy cluster of a given distance matrix.
    
    Parameters
    ----------
    threshold : float
        The threshold that shall be used for the basic clustering of the data.

    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    taxa : list
        An list containing the names of all taxa corresponding to the distances
        in the matrix.

    distances : bool
        If set to ``False``, only the topology of the tree will be returned.

    revert : bool (default=False)
        Specify whether a reverted dictionary should be returned. 

    Notes
    -----
    This is a very simple fuzzy clustering algorithm. It basically does nothing
    else than removing taxa successively from the matrix, flat-clustering the
    remaining taxa with the corresponding threshold, and then returning a
    combined "consensus" cluster in which taxa may be assigned to multiple
    clusters.
    """
    
    g = nx.Graph()
    for taxon in taxa: g.add_node(taxon)

    for idx,taxon in enumerate(taxa):    
    
        new_matrix = []
        for i,line in enumerate(matrix):
            for j,cell in enumerate(line):
                if i < j and i != idx and j != idx:
                    new_matrix += [cell]
        new_matrix = misc.squareform(new_matrix)
        
        clusters = cluster.flat_cluster(
                method,
                threshold,
                new_matrix,
                [t for t in taxa if t != taxon]
                )

        if rcParams['verbose']: print(taxon,idx,clusters)

        for clr in clusters:
            for i,tA in enumerate(clusters[clr]):
                for j,tB in enumerate(clusters[clr]):
                    if i < j:
                        try:
                            g.edge[tA][tB]['weight'] += 1
                        except:
                            g.add_edge(tA,tB,weight=1)
    i = 1
    out = {}
    for c in nx.find_cliques(g):
        out[i] = c
        i += 1

    if revert:
        new_out = {}
        for key,val in out.items():
            for v in val:
                try:
                    new_out[v] += [key]
                except KeyError:
                    new_out[v] = [key]
        return new_out

    return out

def matrix2tree(
        matrix,
        taxa,
        tree_calc = "neighbor",
        distances = True,
        filename = ''
        ):
    """
    Calculate a tree of a given distance matrix.

    Parameters
    ----------
    matrix : {list, numpy.array}
        The distance matrix to be used.
    taxa : list
        A list of the taxa in the distance matrix.
    tree_calc : str (default='neighbor')
        The method for tree calculation that shall be used. Select between:

        * 'neighbor': Neighbor-joining method (:evobib:`Saitou1987`)
        * 'upgma' : UPGMA method (:evobib:`Sokal1958`)
    distances : bool (default=True)
        If set to c{True}, distances will be included in the
        tree-representation.
    filename : str (default='')
        If a filename is specified, the data will be written to that file.

    Returns
    -------
    tree : ~lingpy.thirdparty.cogent.tree.PhyloNode
        A ~lingpy.thirdparty.cogent.tree.PhyloNode object for handling tree
        files.
    """
    
    if tree_calc == 'upgma':
        algorithm = cluster.upgma
    elif tree_calc == 'neighbor':
        algorithm = cluster.neighbor

    newick = algorithm(matrix,taxa,distances)

    tree = cg.LoadTree(treestring=newick)
    
    if not filename:
        return tree
    else:
        out = codecs.open(filename+'.nwk','w','utf-8')
        out.write(str(tree))
        out.close()
        if rcParams['verbose']: print(rcParams['M_file_written'].format(filename,'nwk'))
    
def matrix2groups(
        threshold,
        matrix,
        taxa
        ):
    """
    Calculate flat cluster of distance matrix.

    Parameters
    ----------
    threshold : float
        The threshold to be used for the calculation.
    matrix : {list, numpy.array}
        The distance matrix to be used.
    taxa : list
        A list of the taxa in the distance matrix.

    Returns
    -------
    groups : dict
        A dictionary with the taxa as keys and the group assignment as values.
        
    """
    
    flats = cluster.flat_upgma(
            threshold,
            distances,
            taxa = [t for t in taxa]
            )
    
    mapper = dict(zip(flats,range(1,len(taxa)+1)))
    out = {}
    for key in flats:
        n = 'G_{0}'.format(mapper[key])
        for t in flats[key]:
            out[t] = n
    return out
        

    groups = [flats[i] for i in range(len(taxa))]
    
    # renumber the groups
    groupset = sorted(set(groups))
    renum = dict([(i,j+1) for i,j in zip(
        groupset,
        range(len(groupset))
        )])
    groups = [renum[i] for i in groups]

    return dict(zip(taxa,['G_{0}'.format(g) for g in groups]))

def find_optimal_cutoff(matrix):
    """
    Use the method by :evobib:`Apeltsin2011` in order to find an optimal threshold. 
    """
    
    pass

def link_clustering(
        taxa,
        matrix,
        cutoff=0.5,
        threshold=False,
        revert=False,
        matrix_type = 'distances'
        ):
    """
    Carry out a link clustering analysis using the method by :evobib:`Ahn2010`.

    Parameters
    ----------
    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    taxa : list
        An list containing the names of all taxa corresponding to the distances
        in the matrix.
    
    cutoff : float
        The threshold that shall be used for the initial selection of links
        assignd to the data.


    threshold : float (default=0.5)
        The threshold that shall be used for the internal clustering of the
        data.

    Returns
    -------
    cluster : dict
        A dictionary that displays the clusters.

    """
    # check for matrix type
    if matrix_type == 'distances':
        evaluate = lambda x:True if x < cutoff else False
    elif matrix_type == 'similarities':
        evaluate = lambda x:True if x > cutoff else False
    elif matrix_type == 'weights':
        evaluate = lambda x:False

    # get the edges and the adjacency from the thresholds
    edges = set()
    adjacency = dict([(t,set()) for t in taxa])
    weights = {}

    for i,taxA in enumerate(taxa):
        for j,taxB in enumerate(taxa):
            if i < j:
                if evaluate(matrix[i][j]):
                    edges.add((taxA,taxB))
                    adjacency[taxA].add(taxB)
                    adjacency[taxB].add(taxA)
                elif matrix_type == 'weights':
                    edges.add((taxA,taxB))
                    adjacency[taxA].add(taxB)
                    adjacency[taxB].add(taxA)
                    edges.add((taxB,taxA))
                    weights[taxA,taxB] = matrix[i][j]
                    weights[taxB,taxA] = matrix[i][j]
    
    if not weights:
        weights = None
    
    if edges:
        # initialize the HLC object
        hlc = lc.HLC(adjacency,edges)
    else:
        # check for null edges: if they occur, return the clusters directly
        if revert:
            return dict([(a,[b]) for a,b in zip(taxa,range(len(taxa)))])
        else:
            return dict([(a,[b]) for a,b in zip(range(len(taxa)),taxa)])

    # carry out the analyses using defaults for the clustering
    comms = hlc.single_linkage(threshold=threshold,w=weights)
    edge2cid = comms[0]
    
    # retrieve all clusterings for the nodes
    cluster = dict([(t,[]) for t in taxa])
    
    # retrieve the data
    clr2nodes = dict()
    clr2edges = dict()
    
    for edge,idx in edge2cid.items():
        nodeA,nodeB = edge[0],edge[1]

        try:
            clr2edges[idx] += [edge]
        except KeyError:
            clr2edges[idx] = [edge]
        try:
            clr2nodes[idx] += [nodeA,nodeB]
        except KeyError:
            clr2nodes[idx] = [nodeA,nodeB]


    for idx in clr2nodes:
        clr2nodes[idx] = sorted(set(clr2nodes[idx]))
    
    # delete all clusters that appear as subsets of larger clusters
    delis = []
    for keyA in sorted(clr2nodes):
        for keyB in sorted(clr2nodes):
            if keyA != keyB:
                valsA = set(clr2nodes[keyA])
                valsB = set(clr2nodes[keyB])
                
                if valsA != valsB:
                    if valsA.issubset(valsB):
                        delis += [keyA]
                    elif valsB.issubset(valsA):
                        delis += [keyB]
                elif valsA == valsB:
                    delis += [keyB]
    for k in set(delis):
        del clr2nodes[k]

    # renumber the data
    mapper = dict(zip(clr2nodes.keys(),range(1,len(clr2nodes)+1)))
    
    out = {}
    found = []
    for idx in clr2nodes:
        out[mapper[idx]] = clr2nodes[idx]
        found += clr2nodes[idx]
    missing = [f for f in taxa if f not in found]
    idx = max(out.keys())+1
    for m in missing:
        out[idx] = [m]
        idx += 1

    if not revert:
        return out
    else:
        cluster = dict([(t,[]) for t in taxa])
        for idx in out:
            for t in out[idx]:
                cluster[t] += [idx]

        return cluster

# the following lines of code are devoted to mcl clustering algorithm

def _normalize_matrix(matrix):
    """
    Normalize the matrix.
    """
    return matrix / sum(matrix)

def _is_idempotent(matrix):
    """
    Check whether the matrix is idempotent.
    """

    for line in matrix:
        if len([x for x in set(line) if x != 0]) > 1:
            return False
    return True

def _interprete_matrix(matrix,revert=False):
    """
    Look for attracting nodes in the matrix.
    """

    clusters = []
    flags = len(matrix) * [False]
    for i in range(len(matrix)):
        clr = []
        for j in range(len(matrix)):
            if not flags[j] and matrix[i][j] > 0:
                clr += [j]
                flags[j] = True
        if clr:
            clusters += [clr]
    
    # make a converter for length
    out = [0 for i in range(len(matrix))]
    idx = 1 
    for clr in clusters:
        for i in clr:
            out[i] = idx
        idx += 1
    
    return out

def mcl(
        nodes,
        adjmatrix,
        threshold=False,
        max_steps = 1000,
        inflation = 2,
        expansion = 2,
        add_self_loops = True,
        revert = False
        ):
    """
    Carry out mcl clustering.
    """
    # check for type of matrix
    if type(adjmatrix) != np.ndarray:
        matrix = np.array(adjmatrix)
    else:
        matrix = adjmatrix.copy()

    # check for threshold
    if threshold:
        for i in range(len(matrix)):
            matrix[i][i] = 0
            for j in range(i,len(matrix)):
                if matrix[i][j] > threshold:
                    matrix[i][j] = 0
                    matrix[j][i] = 0
                else:
                    matrix[i][j] = 1
                    matrix[j][i] = 1
    
    # check for self_loops
    if add_self_loops:
        for i in range(len(matrix)):
            matrix[i][i] = sum(matrix[:,i])


    # normalize the matrix
    matrix = _normalize_matrix(matrix)

    # start looping and the like
    steps = 0
    while True:
        
        # expansion
        matrix = np.linalg.matrix_power(matrix,expansion)
        
        # inflation
        matrix = matrix ** inflation

        # normalization
        matrix = _normalize_matrix(matrix)

        # increase steps
        steps += 1

        # check for matrix convergence
        if steps >= max_steps or _is_idempotent(matrix):
            if rcParams['debug']:
                print("[DEBUG] Number of steps {0}.".format(steps))
            break
    
    # retrieve the clusters
    clusters = _interprete_matrix(matrix)

    # modify clusters
    if revert:
        return dict(
                zip(
                    range(len(nodes)),
                    clusters
                    )
                )
    
    return clusters

    
    
        

