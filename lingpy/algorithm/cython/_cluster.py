from __future__ import unicode_literals
try:
    from .misc import transpose,squareform
except ImportError:
    from ._misc import transpose,squareform

def flat_upgma(
        threshold,
        matrix,
        taxa = [],
        revert = False
        ):
    """
    Carry out a flat cluster analysis based on the UPGMA algorithm \
    (:evobib:`Sokal1958`).
    
    Parameters
    ----------
 
    threshold : float
        The threshold which terminates the algorithm.   
    
    matrix : or :py:class:`numpy.array`
        A two-dimensional containing the distances.

    taxa : (default = [])
        A containing the names of the taxa. If the is left empty, the
        indices of the taxa will be returned instead of their names.
    
    Returns
    -------
    
    clusters : dict
        A dictionary with cluster-IDs as keys and a of the taxa
        corresponding to the respective ID as values.

    Examples
    --------
    The function is automatically imported along with LingPy.

    >>> from lingpy import *
    
    Create a of arbitrary taxa.

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
    ~lingpy.algorithm.clustering.upgma
    ~lingpy.algorithm.clustering.neighbor

    """
# [autouncomment]     cdef int i,key
    x = len(matrix)
    out = {}

    clusters = dict([(i,[i]) for i in range(x)])

    tree = []

    _flat_upgma(clusters,matrix,threshold)

    if taxa:
        for key in clusters:
            clusters[key] = [taxa[i] for i in clusters[key]]

        return clusters
    if revert:
        for key in clusters:
            for i in clusters[key]:
                out[i] = key + 1
        return out
    return clusters

def flat_cluster(
        method,
        threshold,
        matrix,
        taxa = [],
        revert = False
        ):
    """
    Carry out a flat cluster analysis based on the UPGMA algorithm.
    
    Parameters
    ----------
    method : { 'upgma', 'single', 'complete' }
        Select between 'ugpma', 'single', and 'complete'.
 
    threshold : float
        The threshold which terminates the algorithm.   
    
    matrix : or :py:class:`numpy.array`
        A two-dimensional containing the distances.

    taxa : (default = [])
        A containing the names of the taxa. If the is left empty, the
        indices of the taxa will be returned instead of their names.
    
    Returns
    -------
    
    clusters : dict
        A dictionary with cluster-IDs as keys and a of the taxa
        corresponding to the respective ID as values.

    Examples
    --------
    The function is automatically imported along with LingPy.

    >>> from lingpy import *
    
    Create a of arbitrary taxa.

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
    ~lingpy.algorithm.clustering.upgma
    ~lingpy.algorithm.clustering.neighbor

    """
# [autouncomment]     cdef int i,key
    x = len(matrix)
    out = {}

    clusters = dict([(i,[i]) for i in range(x)])

    tree = []

    if method == 'upgma':
        _flat_upgma(clusters,matrix,threshold)
    elif method == 'single':
        _flat_single_linkage(clusters,matrix,threshold)
    elif method == 'complete':
        _flat_complete_linkage(clusters,matrix,threshold)

    if taxa:
        for key in clusters:
            clusters[key] = [taxa[i] for i in clusters[key]]

        return clusters
    
    if revert:
        for key in clusters:
            for i in clusters[key]:
                out[i] = key + 1
        return out
    return clusters


def _flat_upgma(
        clusters,
        matrix,
        threshold
        ):
    """
    Internal implementation of flat_upgma.
    """
# [autouncomment]     cdef int i,j,vA,vB,idxA,idxB
# [autouncomment]     cdef list score,valA,valB
    
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
                scores.append(sum(score) / len(score))
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

def _flat_single_linkage(
        clusters,
        matrix,
        threshold
        ):
    """
    Internal implementation of flat_upgma.
    """
# [autouncomment]     cdef int i,j,vA,vB,idxA,idxB
# [autouncomment]     cdef list score,valA,valB
    
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
                scores.append(min(score))
                indices.append((i,j))
    
    minimum = min(scores)
    if minimum <= threshold:
        idxA,idxB = indices[scores.index(minimum)]
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        return _flat_single_linkage(
                clusters,
                matrix,
                threshold
                )
    else:
        pass

def _flat_complete_linkage(
        clusters,
        matrix,
        threshold
        ):
    """
    Internal implementation of flat_upgma.
    """
# [autouncomment]     cdef int i,j,vA,vB,idxA,idxB
# [autouncomment]     cdef list score,valA,valB
    
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
                scores.append(max(score))
                indices.append((i,j))
    
    minimum = min(scores)
    if minimum <= threshold:
        idxA,idxB = indices[scores.index(minimum)]
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        return _flat_complete_linkage(
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

    matrix : or :py:class:`numpy.array`
        A two-dimensional containing the distances.

    taxa : list
        An containing the names of all taxa corresponding to the distances
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
    
    Create an arbitrary of taxa.

    >>> taxa = ['German','Swedish','Icelandic','English','Dutch']
    
    Create an arbitrary matrix.

    >>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])

    Carry out the cluster analysis.

    >>> upgma(matrix,taxa,distances=False)
    '((Swedish,Icelandic),(English,(German,Dutch)));'

    See also
    --------
    ~lingpy.algorithm.clustering.neighbor
    ~lingpy.algorithm.clustering.flat_upgma
   
    """
# [autouncomment]     cdef int i,a,b
    x = len(taxa)
# [autouncomment]     cdef float c,d
# [autouncomment]     cdef str newick_string

    clusters = dict([(i,[i]) for i in range(x)])
    branches = dict([(i,0) for i in range(x)])

    tree = []

    _upgma(clusters,matrix,tree,branches)

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
    
    newick_string = newick[max(newick.keys())] + ';'

    return newick_string

def _upgma(
        clusters,
        matrix,
        tree_matrix,
        branches = {}
        ):
    """
    Internal implementation of the UPGMA algorithm.
    """
# [autouncomment]     cdef int i,vA,vB,idxA,idxB,idxNew
# [autouncomment]     cdef list score,valA,valB

    # check for branches
    if not branches:
        branches = dict([(i,0) for i in clusters])
    
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
                scores.append(sum(score) / len(score))
                indices.append((i,j))

    minimum = min(scores)
    
    idxNew = max(clusters) + 1

    idxA,idxB = indices[scores.index(minimum)]
    
    bA = minimum / 2 - branches[idxA]
    bB = minimum / 2 - branches[idxB]

    branches[idxNew] = minimum / 2

    clusters[idxNew] = clusters[idxA] + clusters[idxB]

    del clusters[idxA]
    del clusters[idxB]

    tree_matrix.append([idxA,idxB,bA,bB])
    
    return _upgma(
            clusters,
            matrix,
            tree_matrix,
            branches
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

    matrix : or :py:class:`numpy.array`
        A two-dimensional containing the distances.

    taxa : list
        An containing the names of all taxa corresponding to the distances
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
    
    Create an arbitrary of taxa.

    >>> taxa = ['Norwegian','Swedish','Icelandic','Dutch','English']
    
    Create an arbitrary matrix.

    >>> matrix = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])

    Carry out the cluster analysis.

    >>> neighbor(matrix,taxa)
    '(((Norwegian,(Swedish,Icelandic)),English),Dutch);'

    See also
    --------
    ~lingpy.algorithm.clustering.upgma
    ~lingpy.algorithm.clustering.flat_upgma
    """
# [autouncomment]     cdef int i,a,b
# [autouncomment]     cdef float c,d
# [autouncomment]     cdef str newick_string
    
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
    
    newick_string = newick[max(newick.keys())] + ';'

    return newick_string

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
# [autouncomment]     cdef int idxA,idxB,idxNew,N,i,j,key
# [autouncomment]     cdef float sAX,sBX,new_score,score,dist_a,dist_b,dist_ab
# [autouncomment]     cdef list line,new_matrix
# [autouncomment]     cdef dict new_clusters
# [autouncomment]     cdef list averages

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
        sAX = matrix[idxA][idxB] / 2
        sBX = matrix[idxA][idxB] / 2

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
        constant_matrix = [[cell for cell in line] for line in matrix]
    
    # get the number of taxa
    N = len(matrix)

    # determine the average scores (divergence r)
    averages = []
    for line in matrix:
        averages.append(sum(line) / (N - 2.0))
    
    # create the new matrix
    new_matrix = [[cell for cell in line] for line in matrix]
    
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
    tmp_score = sum(vals) / len(vals)

    
    # append the indices to the tree matrix
    sAX = matrix[idxA][idxB] / 2.0 + (averages[idxA] - averages[idxB]) / 2
        #(2.0 * (len(matrix) - 2))
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

def _tree2nwk(
        tree,
        taxa,
        distances
        ):
    """
    Convert the tree-matrix created by the _upgma-function to newick representation.
    
    Parameters
    ----------
    tree_matrix : list
        The tree-representation that is yielded by _upgma, also used in
        scipy-cluster algorithms.
    taxa : 
        List of the taxa (or sequences) in the order in which the tree was
        created.
    distances : bool
        Specify whether distances should be included in the string or not.

    Returns
    -------
    newick : str
        A newick string.

    """

# [autouncomment]     cdef int i,a,b
# [autouncomment]     cdef float c,d
    x = len(taxa)
# [autouncomment]     cdef str newick_string

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
    
    newick_string = newick[max(newick.keys())] + ';'

    return newick_string    
    
    
    
