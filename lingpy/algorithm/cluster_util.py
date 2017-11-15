"""Various utility functions which are useful for algorithmic operations"""

from itertools import combinations_with_replacement
from collections import defaultdict
import random

def valid_cluster(sequence):
    """Only allow to have sequences which have consecutive ordering of elements.

    Parameters
    ----------
    sequence : list
        A cluster sequence in which elements should be consecutively ordered, starting from
        0, and identical segments in the sequence retrieve the same number.

    Returns
    -------
    valid_cluster : bool
        True, if the cluster is valid, and False if it judged to be invalid.

    Examples
    --------
    We define a valid and an invalid cluster sequence:

        >>> clrA = [0, 1, 2, 3]
        >>> clrB = [1, 1, 2, 3] # should be [0, 0, 1, 2]
        >>> from lingpy.algorithm.utils import valid_cluster
        >>> valid_cluster(clrA)
        True
        >>> valid_cluster(clrB)
        False

    Seealso:
    --------
    generate_all_clusters
    generate_random_cluster
    order_cluster
    mutate_cluster

    """
    visited = set()
    maxelm = -1
    for segment in sequence:
        if segment not in visited:
            if segment != maxelm + 1:
                return False
            visited.add(segment)
            maxelm = segment
    return True


def generate_all_clusters(numbers):
    """Generate all possible clusters for a number of elements.
    
    Returns
    -------
    clr : iterator
        An iterator that will yield the next of all possible clusters.

    Seealso
    -------
    valid_cluster
    generate_random_cluster
    order_cluster
    mutate_cluster
    """

    for clr in combinations_with_replacement(range(numbers), numbers):
        if valid_cluster(clr):
            yield clr


def generate_random_cluster(numbers, bias=False):
    """Generate a random cluster for a number of elements.
    
    Parameters
    ----------
    numbers : int
        Number of separate entities which should be clustered.
    bias : str (default=False)
        When set to "lumper" will tend to create larger groups, when set to
        "splitter" it will tend to produce smaller groups.

    Returns
    -------
    cluster : list
        A list with consecutive ordering of clusters, starting from zero.

    Seealso
    -------
    valid_cluster
    generate_all_clusters
    order_cluster
    mutate_cluster

    """
    out = []
    maxelm = 0
    if bias == 'lumper':
        selector = lambda x, y, z: z * 2 + (abs(x-y)+1) * [x]
    elif bias == 'splitter':
        selector = lambda x, y, z: ((x-y)+1) * [x]

    for i in range(numbers):
        if bias:
            select_from = []
            for i in range(maxelm+1):
                select_from += selector(i, maxelm, [o for o in out])
            nextelm = random.choice(select_from)
        else:
            nextelm = random.randint(0, maxelm)
        if nextelm == maxelm:
            maxelm += 1
        out += [nextelm]
    return out


def order_cluster(clr):
    """Order a cluster into the form of a valid cluster.
    
    Parameters
    ----------
    clr : list
        A list with clusters assigned by given each element a specific clusuter
        ID.

    Returns
    -------
    valid_cluster : list
        A list in which the IDs start from zero and increase consecutively with
        each new cluster introduced.

    Seealso
    -------
    valid_cluster
    generate_all_clusters
    generate_random_cluster
    mutate_cluster
    """
    current = 0
    dct = {}
    out = []
    for element in clr:
        pendant = dct.get(element, current)
        out += [pendant]
        if element not in dct:
            dct[element] = current
            current += 1
    return out


def mutate_cluster(clr, chance=0.5):
    """Mutate a cluster.
    
    Parameters
    ----------
    clr : cluster
        A list with ordered clusters.
    chance : float (default=0.5)
        The mutation rate for each element in a cluster. If set to 0.5, this
        means that in 50% of the cases, an element will be assigned to another
        cluster or a new cluster.

    Returns
    -------
    valid_cluster : list
        A newly clustered list in consecutive order.
    
    Seealso
    -------
    valid_cluster
    generate_all_clusters
    generate_random_cluster
    order_cluster

    """
    out = []
    numbers = list(range(len(clr)))
    for element in clr:
        if random.uniform(0, 1) < chance:
            out += [random.choice(numbers)]
        else:
            out += [element]

    return order_cluster(out)



