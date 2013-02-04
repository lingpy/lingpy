"""
This module provides functions to caculate basic distance measures.

"""
__author__ = "Steven Moran"
__date__ = "2013-2"

import itertools
import operator

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
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(map(ne, str1, str2))

def jaccard(set1, set2):
    """
    Computer the Jaccard distance between two sets.

    Jaccard distance measures the dissimilarity between sample sets.
    It is complementary to the Jaccard coefficient and is obtained by 
    subtracting the Jaccard coefficient from 1, or, equivalently, by 
    dividing the difference of the sizes of the union and the intersection 
    of two sets by the size of the union:

    J(A,B) = 1 - J(A,B) = |A ∪ B| - |A ∩ B| / |A ∪ B|

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
    assert(type(set1) and type(set2) == set)
    n = len(set1.intersection(set2))
    return (n/float(len(set1) + len(set2) -n))

def euclidean():
    pass

def pearson():
    pass

if __name__=="__main__":
    print("test Hamming:")
    str1 = "aaa"
    str2 = "aba"
    print("str1 = ", str1)
    print("str2 = ", str2)
    print("hamming distances = ", hamming(str1, str2))
    print()

    print("test Jaccard:")
    set1 = set("abcd")
    set2 = set("cdef")
    print("set1 = ", set1)
    print("set2 = ", set2)
    print("jaccard distance = ", jaccard(set1, set2))
    print()
