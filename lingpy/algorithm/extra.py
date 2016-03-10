"""
Adapting specific cluster algorithms from scikit-learn to LingPy.
"""
from lingpy import log

try:
    from sklearn import cluster
except ImportError:
    log.missing_module('sklearn')
import numpy as np


def dbscan(
        threshold,
        matrix,
        taxa,
        revert=False,
        min_samples=1):
    """
    Compute DBSCAN cluster analysis.
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

    clr = {}
    for i, t in enumerate(taxa):
        try:
            clr[labels[i]] += [t]
        except KeyError:
            # FIXME: clusters is undefined!
            clr[clusters[i]] = [t]

    return clr


def affinity_propagation(threshold, matrix, taxa, revert=False):
    """
    Compute affinity propagation from the matrix.
    """
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
                matrix[i][j] = - np.log2(1 - score ** 2)  # -np.log2(score+0.01)
                matrix[j][i] = matrix[i][j]  # score ** 2  # -np.log2(score+0.01)
            else:
                matrix[i][j] = - score ** 5  # 0.0
                matrix[j][i] = - score ** 5  # 0.0

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

    clr = {}
    for i, t in enumerate(taxa):
        try:
            clr[labels[i]] += [t]
        except KeyError:
            # FIXME: clusters is undefined!
            clr[clusters[i]] = [t]

    return clr
