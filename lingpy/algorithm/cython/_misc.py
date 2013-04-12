# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-13 08:18
# modified : 2013-03-13 08:18
"""
Module provides miscellaneous functions and routines.
"""

__author__="Johann-Mattis List"
__date__="2013-03-13"


#cdef extern from "math.h": 
from numpy import sqrt

def transpose(
        matrix
        ):
    """
    Transpose a matrix along its two dimensions.

    Parameters
    ----------
    matrix : list
        A two-dimensional list.
    """
#     cdef int i,j
    lA = len(matrix)
    lB = len(matrix[0])

    out = [[matrix[i][j] for i in range(lA)] for j in range(lB)]

    return out

def squareform(
        x
        ):
    """
    A simplified version of the :py:func:`scipy.spatial.distance.squareform` \
    function.

    Parameters
    ----------

    x : :py:class:`numpy.array` or list
        The one-dimensional flat representation of a symmetrix distance matrix.

    Returns
    -------
    matrix : :py:class:`numpy.array`
        The two-dimensional redundant representation of a symmetric distance matrix.

    """
#     cdef int i,j,k
    l = len(x)

    # calculate the length of the square
    s = int(sqrt(2 * l) + 1)
    
    out = [[0.0 for i in range(s)] for j in range(s)]
    
    k = 0
    for i in range(s):
        for j in range(s):
            if i < j:
                out[i][j] = x[k]
                out[j][i] = x[k]
                k += 1
    return out

class ScoreDict(object):
    """
    Class allows quick access to scoring functions with using dictionary \
    syntax.

    Parameters
    ----------
    chars : list
        The list of all character tokens for the scoring dictionary.
    matrix : list
        A two-dimensional scoring matrix.

    """
    def __init__(
            self,
            chars,
            matrix
            ):
#         cdef int i
#         cdef str character

        self.chars2int = dict([(character,i) for character,i in
            zip(chars,range(len(chars)))])

        self.matrix = matrix

    def __getitem__(
            self,
            x
            ):
#         cdef int i,j
        try:
            i = self.chars2int[x[0]]
            j = self.chars2int[x[1]]

            return self.matrix[i][j]
        except:
            return -22.5

    def __setitem__(
            self,
            x,
            y
            ):
#         cdef int i,j
        i = self.chars2int[x[0]]
        j = self.chars2int[x[1]]

        if i == j:
            self.matrix[i][j] = y
            self.matrix[j][i] = y
        else:
            self.matrix[i][i] = y

    def __repr__(self):
        return str(list(self.chars2int.keys()))

    def __str__(self):
        return str(list(self.chars2int.items()))

