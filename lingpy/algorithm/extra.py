"""
Adapting specific cluster algorithms from scikit-learn to LingPy.
"""
from collections import defaultdict
from lingpy import log


try:
    from sklearn import cluster
except ImportError:
    log.missing_module('sklearn')
    cluster = False
try:
    import igraph
except ImportError:
    log.missing_module('igraph')
    igraph = False

import numpy as np

def dbscan(
        threshold,
        matrix,
        taxa,
        revert=False,
        min_samples=1):
    """
    Compute DBSCAN cluster analysis.

    Parameters
    ----------
    threshold : float
        The threshold for clustering you want to use.
    matrix : list
        The two-dimensional matrix passed as list or array.
    taxa : list
        The list of taxon names. If set to "False" a fake list of taxon names
        will be created, giving a positive numerical ID in increasing order for
        each column in the matrix.
    revert : bool
        If set to "False", don't return taxon names but simply the language
        identifiers and their labels as a dictionary. Otherwise returns a
        dictionary with labels as keys and list of taxon names as values.
    min_samples : int (default=1)
        The minimal samples parameter of the DBCSCAN method from the SKLEARN
        package.
    
    Returns
    -------
    clusters : dict
        Either a dictionary of taxon identifiers and labels, or a dictionary of
        labels and taxon names.

    Notes
    -----
    This method does not work as expected, probably since it normally requires
    distances between points as input. We list it only for completeness here,
    but urge to be careful when using the code and checking properly our
    implementation in the source code.
    
    Requires the scikitlearn package, downloadable from http://scikit-learn.org/.
    """
    if not taxa:
        taxa = list(range(1, len(matrix) + 1))

    core_samples, labels = cluster.dbscan(
        matrix, eps=threshold, min_samples=min_samples, metric='precomputed')
    
    # change to our internal cluster style
    idx = max(labels) + 1
    if idx == 0:
        idx += 1
    for i, c in enumerate(labels):
        if c == -1:
            labels[i] = idx
            idx += 1

    # check for revert
    if revert:
        return dict(zip(range(len(taxa)), labels))

    clr = defaultdict(list)
    for i, t in enumerate(taxa):
        clr[labels[i]] += [t]
    return clr

def affinity_propagation(threshold, matrix, taxa, revert=False):
    """
    Compute affinity propagation from the matrix.
    
    Parameters
    ----------
    threshold : float
        The threshold for clustering you want to use.
    matrix : list
        The two-dimensional matrix passed as list or array.
    taxa : list
        The list of taxon names. If set to "False" a fake list of taxon names
        will be created, giving a positive numerical ID in increasing order for
        each column in the matrix.
    revert : bool
        If set to "False", don't return taxon names but simply the language
        identifiers and their labels as a dictionary. Otherwise returns a
        dictionary with labels as keys and list of taxon names as values.

    Returns
    -------
    clusters : dict
        Either a dictionary of taxon identifiers and labels, or a dictionary of
        labels and taxon names.
    
    Notes
    -----

    Affinity propagation is a clustering method originally proposed by
    :evobib:`Frey2007`.

    Requires the scikitlearn package, downloadable from http://scikit-learn.org/.



    """
    if not cluster:
        raise ValueError("The package scikitlearn is needed to run this analysis.")
    if not taxa:
        taxa = list(range(1, len(matrix) + 1))
    # turn distances to similarities
    matrix = np.array(matrix)

    # iterate over matrix
    for i, line in enumerate(matrix):
        matrix[i][i] = 10
        for j in range(i + 1, len(matrix)):
            score = matrix[i][j]
            if score < threshold:
                matrix[i][j] = - np.log2(1 - score ** 2)  
                matrix[j][i] = matrix[i][j]  
            else:
                matrix[i][j] = - score ** 5
                matrix[j][i] = - score ** 5

    ap = cluster.AffinityPropagation(affinity='precomputed')
    labels = ap.fit_predict(matrix)

    # change to our internal cluster style
    idx = max(labels) + 1
    if idx == 0:
        idx += 1
    for i, c in enumerate(labels):
        if c == -1:
            labels[i] = idx
            idx += 1

    # check for revert
    if revert:
        return dict(zip(range(len(taxa)), labels))

    clr = defaultdict(list)
    for i, t in enumerate(taxa):
        clr[labels[i]] += [t]
    return clr

def infomap_clustering(threshold, matrix, taxa=False, revert=False):
    """
    Compute the Infomap clustering analysis of the data.

    Parameters
    ----------
    threshold : float
        The threshold for clustering you want to use.
    matrix : list
        The two-dimensional matrix passed as list or array.
    taxa : list
        The list of taxon names. If set to "False" a fake list of taxon names
        will be created, giving a positive numerical ID in increasing order for
        each column in the matrix.
    revert : bool
        If set to "False", don't return taxon names but simply the language
        identifiers and their labels as a dictionary. Otherwise returns a
        dictionary with labels as keys and list of taxon names as values.

    Returns
    -------
    clusters : dict
        Either a dictionary of taxon identifiers and labels, or a dictionary of
        labels and taxon names.

    Notes
    -----
    Infomap clustering is a community detection method originally proposed by
    :evobib:`Rosvall2008`.

    Requires the igraph package is required, downloadable from http://igraph.org/.
    """
    if not igraph:
        raise ValueError("The package igraph is needed to run this analysis.")
    if not taxa:
        taxa = list(range(1, len(matrix) + 1))

    G = igraph.Graph()
    vertex_weights = []
    for i in range(len(matrix)):
        G.add_vertex(i)
        vertex_weights += [0]
    
    # variable stores edge weights, if they are not there, the network is
    # already separated by the threshold
    for i,row in enumerate(matrix):
        for j,cell in enumerate(row):
            if i < j:
                if cell <= threshold:
                    G.add_edge(i, j)
     
    comps = G.community_infomap(edge_weights=None,
            vertex_weights=None)
    D = {}
    for i,comp in enumerate(comps.subgraphs()):
        vertices = [v['name'] for v in comp.vs]
        for vertex in vertices:
            D[vertex] = i+1

    if revert:
        return D

    clr = defaultdict(list)
    for i,t in enumerate(taxa):
        clr[D[i]] += [t]
    return clr
