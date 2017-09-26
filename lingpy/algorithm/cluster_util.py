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
    order_cluster

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
    """Generate all possible clusters for a number of elements."""
    for clr in combinations_with_replacement(range(numbers), numbers):
        if valid_cluster(clr):
            yield clr


def generate_random_cluster(numbers):
    """Generate a random cluster for a number of elements."""
    out = []
    maxelm = 0
    for i in range(numbers):
        nextelm = random.randint(0,maxelm)
        if nextelm == maxelm:
            maxelm += 1
        out += [nextelm]
    return out


def order_cluster(clr):
    """Order a cluster into our numeric form."""
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
    """Mutate a cluster."""
    out = []
    numbers = list(range(len(clr)))
    for element in clr:
        if random.uniform(0, 1) < chance:
            out += [random.choice(numbers)]
        else:
            out += [element]

    return order_cluster(out)



