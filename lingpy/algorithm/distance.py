# *-* coding: utf-8 *-*
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
"""
This module provides functions to calculate basic distance measures.

"""
__author__ = "Steven Moran"
__date__ = "2013-2"

import operator

from six.moves import zip_longest


def hamming(str1, str2):
    """
    Compute the Hamming distance between two strings.

    The Hamming distance (see :evobib:`Hamming1950`) is defined as the
    number of bits that are different between two vectors.

    Parameters
    ----------
    str1 : str
        str to be compared to str2

    str2 : str
        str to be compared to str1

    Returns
    -------
    _ : int
        the hamming distance 

    """
    return sum(operator.ne(*pair) for pair in zip_longest(str1, str2, fillvalue=None))


def jaccard(set1, set2):
    """
    Computer the Jaccard distance between two sets.

    Jaccard distance measures the dissimilarity between sample sets.
    It is complementary to the Jaccard coefficient and is obtained by 
    subtracting the Jaccard coefficient from 1, or, equivalently, by 
    dividing the difference of the sizes of the union and the intersection 
    of two sets by the size of the union:

    J(A,B) = 1 - J(A,B) = \|A ∪ B\| - \|A ∩ B\| / \|A ∪ B\|

    Parameters
    ----------
    set1 : set
        set to be compared to set2

    set2 : set
        set to be compared to set1

    Returns
    -------
    _ : float
        the Jaccard distance

    """
    assert isinstance(set1, set) and isinstance(set2, set)
    if not set1 and not set2:
        return 0
    n = len(set1.intersection(set2))
    return n / float(len(set1) + len(set2) - n)


def euclidean(p, q):
    """
    Distance between two points (p, q) in any dimension of space.

    """
    sum_of_squares = 0.0

    for i in range(len(p)):
        sum_of_squares += (p[i] - q[i]) ** 2

    return sum_of_squares ** 0.5


def pearson(x, y):
    """
    Pearson correlation coefficient.

    n = len(x)
    values = range(n)

    sum_x = sum([float(x[i]) for i in values])
    sum_y = sum([float(y[i]) for i in values])

    sqr_sum_x = sum([x[i]**2.0 for i in values])
    sqr_sum_y = sum([y[i]**2.0 for i in values])

    p_sum = ([x[i]*y[i] for i in values])

    num = p_sum-(sum_x * sum_y / n)
    den = ((sqr_sum_x - pow(sum_x,2)/n) * (sqr_sum_y - pow(sum_y,2)/n))**.5

    if den == 0:
        return 0

    r = dum/den

    return r
    """
    raise NotImplementedError()
