cdef extern from "math.h": 
    double sqrt( double x)

def transpose(
        list matrix
        ):
    """
    Transpose a matrix along its two dimensions.

    Parameters
    ----------
    matrix : list
        A two-dimensional list.
    """
    cdef int i,j
    cdef int lA = len(matrix)
    cdef int lB = len(matrix[0])

    cdef list out = [[matrix[i][j] for i in range(lA)] for j in range(lB)]

    return out

def squareform(
        list x
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
    cdef int i,j,k
    cdef int l = len(x)

    # calculate the length of the square
    cdef int s = int(sqrt(2 * l) + 1)
    
    cdef list out = [[0.0 for i in range(s)] for j in range(s)]
    
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
    Class allows quick access to scoring functions using dictionary \
    syntax.

    Parameters
    ----------
    chars : list
        The list of all character tokens for the scoring dictionary.
    matrix : list
        A two-dimensional scoring matrix.

    Notes
    -----
    Since this class has dictionary syntax, you can always also just create a
    dictionary in order to store your scoring functions. Scoring dictionaries
    should contain a tuple of segments to be compared as a key, and a float or
    integer as a value, with negative values indicating dissimilarity, and
    positive values similarity.

    Examples
    --------
    Initialize a ScoreDict object::
        >>> from lingpy.algorith.cython.misc import ScoreDict
        >>> scorer = ScoreDict(['a', 'b'], [1, -1, -1, 1])
        
    Retrieve scores::
        >>> scorer['a', 'b']
        -1
        >>> scorer['a', 'a']
        1
        >>> scorer['a', 'X']
        -22.5

    """
    def __init__(
            self,
            list chars,
            list matrix
            ):
        cdef int i
        cdef str character

        self.chars2int = dict([(character,i) for character,i in
            zip(chars,range(len(chars)))])

        self.matrix = matrix

    def __getitem__(
            self,
            tuple x
            ):
        cdef int i,j
        try:
            i = self.chars2int[x[0]]
            j = self.chars2int[x[1]]

            return self.matrix[i][j]
        except:
            return -22.5

    def __setitem__(
            self,
            tuple x,
            float y
            ):
        cdef int i,j
        i = self.chars2int[x[0]]
        j = self.chars2int[x[1]]

        if i != j:
            self.matrix[i][j] = y
            self.matrix[j][i] = y
        else:
            self.matrix[i][i] = y

    def __repr__(self):
        return str(list(self.chars2int.keys()))

    def __str__(self):
        return str(list(self.chars2int.items()))

