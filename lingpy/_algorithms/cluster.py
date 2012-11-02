"""
This module provides functions for basic cluster algorithms.

"""
from __future__ import division,print_function
from numpy import array
from misc import squareform


def flat_upgma(
        matrix,
        threshold,
        taxa = None
        ):
    """
    Carry out a flat cluster analysis based on the UPGMA algorithm \
    (:evobib:`Sokal1958`).
    
    Parameters
    ----------
    
    matrix : list or :py:class:`numpy.array`
        A two-dimensional list containing the distances.

    threshold : float
        The threshold which terminates the algorithm.

    taxa : list
        A list containing the names of the taxa. If set to ``None``, the
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

    >>> flat_upgma(clusters,matrix,0.5)
    {0: ['German', 'Dutch', 'English'], 1: ['Swedish', 'Icelandic']}

    See also
    --------
    lingpy.algorithm.clusters.upgma
    lingpy.algorithm.clusters.neighbor

    """
    x = len(taxa)

    clusters = dict([(i,[i]) for i in range(x)])

    tree = []

    _flat_upgma(clusters,matrix,threshold)

    if taxa:
        for key in clusters:
            clusters[key] = [taxa[i] for i in clusters[key]]

        return clusters
    
    return clusters

def _flat_upgma(
        clusters,
        matrix,
        threshold
        ):
    """
    Internal implementation of flat_upgma.
    """
    
    # terminate when the dictionary is of length 1
    if len(clusters) == 1:
        return

    scores = []
    indices = []

    for i,valA in clusters.items():
        for j,valB in clusters.items():
            if i != j:
                score = []
                for vA in valA:
                    for vB in valB:
                        score += [matrix[vA][vB]]
                score = array(score)
                score = score.mean() - 0.25 * score.std()
                scores.append(score)
                indices.append((i,j))

    minimum = min(scores)
    if minimum <= threshold:
        idxA,idxB = indices[scores.index(minimum)]
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        return _flat_upgma(
                clusters,
                matrix,
                threshold
                )
    else:
        pass

def upgma(
        matrix,
        taxa,
        distances = True,
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
    x = len(taxa)

    clusters = dict([(i,[i]) for i in range(x)])

    tree = []

    _upgma(clusters,matrix,tree)

    newick = dict([(i,taxa[i]) for i in range(len(taxa))])
    
    # create different output, depending on the options for the inclusion of
    # distances or topology only
    if distances:
        for i,(a,b,c,d) in enumerate(tree):
            newick[x+i] = '({0}:{2:.2f},{1}:{3:.2f})'.format(
                    newick[a],
                    newick[b],
                    c,
                    d
                    )
    else:
        for i in range(len(tree)):
            newick[x+i] = '({0},{1})'.format(
                    newick[tree[i][0]],
                    newick[tree[i][1]]
                    )
    
    newick = newick[max(newick.keys())] + ';'

    return newick

def _upgma(
        clusters,
        matrix,
        tree_matrix
        ):
    """
    Internal implementation of the UPGMA algorithm.
    """
    
    # terminate when the dictionary is of length 1
    if len(clusters) == 1:
        return

    scores = []
    indices = []

    for i,valA in clusters.items():
        for j,valB in clusters.items():
            if i != j:
                score = []
                for vA in valA:
                    for vB in valB:
                        score += [matrix[vA][vB]]
                score = array(score)
                score = score.mean()
                scores.append(score)
                indices.append((i,j))

    minimum = min(scores)
    
    idxNew = max(clusters) + 1

    idxA,idxB = indices[scores.index(minimum)]
    
    clusters[idxNew] = clusters[idxA] + clusters[idxB]

    del clusters[idxA]
    del clusters[idxB]

    tree_matrix.append([idxA,idxB,minimum/2,minimum/2])
    
    return _upgma(
            clusters,
            matrix,
            tree_matrix
            )

def neighbor(
        matrix,
        taxa,
        distances = True
        ):
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
    
    x = len(taxa)

    clusters = dict([(i,[i]) for i in range(x)])
    tree = []

    _neighbor(clusters,matrix,tree)

    newick = dict([(i,taxa[i]) for i in range(x)])
        
    # create different output, depending on the options for the inclusion of
    # distances or topology only
    if distances:
        for i,(a,b,c,d) in enumerate(tree):
            newick[x+i] = '({0}:{2:.2f},{1}:{3:.2f})'.format(
                    newick[a],
                    newick[b],
                    c,
                    d
                    )
    else:
        for i in range(len(tree)):
            newick[x+i] = '({0},{1})'.format(
                    newick[tree[i][0]],
                    newick[tree[i][1]]
                    )
    
    newick = newick[max(newick.keys())] + ';'

    return newick

def _neighbor(
        clusters,
        matrix,
        tree_matrix,
        constant_matrix = [],
        tracer = {}
        ):
    """
    Internal implementation of the neighbor-joining algorithm.
    """
    if len(clusters) == 1:
        return

    # define a tracer for the order of the tree
    if not tracer:
        tracer = dict([(tuple([a]),b[0]) for (a,b) in clusters.items()])
    # terminate when the dictionary is of length 1
    if len(clusters) == 2:
        
        idxA,idxB = 0,1

        # create the new index for the tracer
        idxNew = max(tracer.values()) + 1

        tracer[tuple(clusters[idxA]+clusters[idxB])] = idxNew

        # append the indices to the tree matrix
        sAX = matrix[idxA][idxB] 
        sBX = matrix[idxA][idxB]

        tree_matrix.append(
                (
                    tracer[tuple(clusters[idxA])],
                    tracer[tuple(clusters[idxB])],
                    sAX,
                    sBX
                    )
                )

        # join the clusters according to the index
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        
        return

    # create the constant matrix when the process starts
    if not constant_matrix:
        constant_matrix = list(matrix.copy())
    
    # get the number of taxa
    N = len(matrix)

    # determine the average scores (divergence r)
    averages = []
    for line in matrix:
        averages.append(sum(line) / (N - 2.0))
    
    # create the new matrix
    new_matrix = matrix.copy()
    
    # fill in the new scores
    for i,line in enumerate(matrix):
        for j,score in enumerate(line):
            if i > j:
                new_score = score - averages[i] - averages[j]
                new_matrix[i][j] = new_score
                new_matrix[j][i] = new_score
    
    # determine the minimal score
    scores = []
    indices = []
    for i in sorted(clusters.keys()):
        for j in sorted(clusters.keys()):
            if i < j:
                scores.append(new_matrix[i][j])
                indices.append((i,j))

    minimum = min(scores)
    idxA,idxB = indices[scores.index(minimum)]

    # check for the average of the clusters
    vals = []
    for i in clusters[idxA]:
        for j in clusters[idxB]:
            vals.append(constant_matrix[i][j])
    tmp_score = array(vals).mean()

    
    # append the indices to the tree matrix
    sAX = (matrix[idxA][idxB] + averages[idxA] - averages[idxB]) / 2.0
    sBX = matrix[idxA][idxB] - sAX
    tree_matrix.append(
            (
                tracer[tuple(clusters[idxA])],
                tracer[tuple(clusters[idxB])],
                sAX,
                sBX
                )
            )

    # create the new index for the tracer
    idxNew = max(tracer.values()) + 1
    tracer[tuple(clusters[idxA]+clusters[idxB])] = idxNew

    # join the clusters according to the index
    clusters[idxA] += clusters[idxB]
    del clusters[idxB]

    # create new cluster-dictionary
    new_clusters = {}

    # insert values in new dictionary
    for i,key in enumerate(sorted(clusters.keys())):
        new_clusters[i] = clusters[key]

    # create new matrix
    new_matrix = []

    # iterate over old matrix and fill in keys for new matrix
    dist_ab = matrix[idxA][idxB]
    for i,a in enumerate(sorted(clusters.keys())):
        for j,b in enumerate(sorted(clusters.keys())):
            if i < j and a != idxA and b != idxA:
                new_matrix.append(matrix[a][b])
            elif i < j and a == idxA:
                dist_a = matrix[idxA][b]
                dist_b = matrix[idxB][b]
                new_matrix.append(((dist_a + dist_b) - dist_ab) / 2.0)
            elif i < j and b == idxA:
                dist_a = matrix[idxA][a]
                dist_b = matrix[idxB][a]
                new_matrix.append(((dist_a + dist_b) - dist_ab) / 2.0)
    
    # get values of new_clusters into clusters
    clusters = {}
    for key,val in new_clusters.items():
        clusters[key] = val

    # return score
    return _neighbor(
            clusters,
            squareform(new_matrix),
            tree_matrix,
            constant_matrix,
            tracer
            )

def _flat_neighbor(
        matrix,
        taxa,
        threshold
        ):
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
    
    x = len(taxa)

    clusters = dict([(i,[i]) for i in range(x)])
    tree = []

    _upgma(clusters,matrix,tree)

    newick = dict([(i,[i]) for i in range(x)])
    newick[x] = True
    clusters = {}
    y = 0
    # create different output, depending on the options for the inclusion of
    # distances or topology only
    for i,(a,b,c,d) in enumerate(tree):
        if newick[a] and newick[b]:
            scores = []
            for A in newick[a]:
                for B in newick[b]:
                    scores.append(matrix[A][B])
            score = array(scores)
            score = score.mean() - 0.25 * score.std()

            if score <= threshold and i != len(tree) -1:
                newick[x+i] = newick[a] + newick[b]
            elif score > threshold and i != len(tree) - 1:
                for v in newick[a]:
                   clusters[v] = y
                y += 1
                for v in newick[b]:
                    clusters[v] = y
                y += 1
                newick[x+i] = False
            elif score <= threshold and i == len(tree) - 1:
                for v in newick[a] + newick[b]:
                    clusters[v] = y
            elif score > threshold and i == len(tree) - 1:
                for v in newick[a]:
                   clusters[v] = y
                y += 1
                for v in newick[b]:
                    clusters[v] = y
                y += 1
        else:
            newick[x+i] = False
            if newick[a]:
                for v in newick[a]:
                    clusters[v] = y
                y += 1
            if newick[b]:
                for v in newick[b]:
                    clusters[v] = y
                y += 1
    if clusters:
        out = {}
        for i,j in clusters.items():
            try:
                out[j] += [taxa[i]]
            except KeyError:
                out[j] = [taxa[i]]
    else:
        out = dict([(i,[taxa[i]]) for i in range(x)])
    
    return out
    

