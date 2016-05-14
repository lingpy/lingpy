# we start with basic alignment functions
def globalign(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer
        ):
    """
    Carry out global alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    
    Notes
    -----
    This is the function that is called to carry out global alignment analyses
    when using many of LingPy's classes for alignment analyses, like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. It differs from classical
    Needleman-Wunsch alignment (compare :evobib:`Needleman1970`) in a couple of aspects.
    These include, among others, the use of a gap extension *scale* rather than
    a gap extension penalty (the scale consecutively reduces the gap penalty
    and thus lets gap penalties approach zero if gapped regions are large), the
    use of individual gap opening penalties for all positions of a sequence,
    and the use of prosodic strings, and prosodic factors that raise scores
    when segments occur in the same prosodic environment. 

    If one sets certain of these parameters to zero or one and uses the same
    gap opening penalties, however, the function will
    behave like the traditional Needleman-Wunsch algorithm, and since it is
    implemented in Cython, it will work faster than a pure Python
    implementation for alignment algorithms.
    
    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score.

    Examples
    --------
    We show that the Needleman-Wunsch algorithms yields the same result as the
    globalign algorithm, provided we adjust the parameters::
        >>> from lingpy.algorithm.cython.calign import globalign
        >>> from lingpy.align.pairwise import nw_align
        >>> nw_align('abab', 'baba')
        (['a', 'b', 'a', 'b', '-'], ['-', 'b', 'a', 'b', 'a'], 1)

        >>> globalign(list('abab'), list('baba'), 4 * [-1], 4 * [-1], 'aaaa', 'aaaa', 4, 4, 1, 0, {("a","b"):-1, ("b","a"): -1, ("a","a"): 1, ("b", "b"): 1})
        (['a', 'b', 'a', 'b', '-'], ['-', 'b', 'a', 'b', 'a'], 1.0)
    
    See also
    --------
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.dialign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef float gapA, gapB, match, sim

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # modify matrix and traceback
    traceback[0][0] = 1
    for i in range(1,M+1):
        matrix[0][i] = matrix[0][i-1] + gopA[i-1] * scale
        traceback[0][i] = 2
    for i in range(1,N+1):
        matrix[i][0] = matrix[i-1][0] + gopB[i-1] * scale
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]

            # calculate costs for gapB
            if traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]

            # calculate costs for match

            # get the score
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match += matrix[i-1][j-1] + match * factor
            elif abs(ord(proA[j-1])-ord(proB[i-1])) >= 2:
                match += matrix[i-1][j-1] + match * factor / 2
            else:
                match += matrix[i-1][j-1]

            # determine minimal cost
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[N][M]

    # carry out the traceback
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA += ['-']
            almB += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA += [seqA[j-1]]
            almB += [seqB[i-1]]
            i -= 1
            j -= 1
        else:
            almA += [seqA[j-1]]
            almB += ['-']
            j -= 1

    # turn alignments back
    almA,almB = almA[::-1],almB[::-1]

    # return alignments
    return almA,almB,sim

def secondary_globalign(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer,
        str r # restricted_chars
        ):
    """
    Carry out global alignment of two sequences with secondary sequence structures.

    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    r : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
    
    Notes
    -----
    This is the function that is called to carry out global alignment analyses
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. It differs from classical
    Needleman-Wunsch alignment (compare :evobib:`Needleman1970`) in a couple of aspects.
    These include, among others, the use of a gap extension *scale* rather than
    a gap extension penalty (the scale consecutively reduces the gap penalty
    and thus lets gap penalties approach zero if gapped regions are large), the
    use of individual gap opening penalties for all positions of a sequence,
    and the use of prosodic strings, and prosodic factors that raise scores
    when segments occur in the same prosodic environment. 

    If one sets certain of these parameters to zero or one and uses the same
    gap opening penalties, however, the function will
    behave like the traditional Needleman-Wunsch algorithm, and since it is
    implemented in Cython, it will work faster than a pure Python
    implementation for alignment algorithms.

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score.

    Examples
    --------
    We compare globalign with secondary_globalign::
        >>> from lingpy.algorithm.cython.calign import globalign, secondary_globalign
        >>> globalign(list('abab'), list('baba'), 4 * [-1], 4 * [-1], 'aaaa', 'aaaa', 4, 4, 1, 0, {("a","b"):-1, ("b","a"): -1, ("a","a"): 1, ("b", "b"): 1})
        (['a', 'b', 'a', 'b', '-'], ['-', 'b', 'a', 'b', 'a'], 1.0)
        >>> secondary_globalign(list('ab.ab'), list('ba.ba'), 5 * [-1], 5 * [-1], 'ab.ab', 'ba.ba', 5, 5, 1, 0, {("a","b"):-1, ("b","a"): -1, ("a","a"): 1, ("b", "b"): 1, ("a",".") : -1, ("b","."):-1, (".","."):0, (".", "b"): -1, (".", "a"):-1}, '.')
        (['a', 'b', '-', '.', 'a', 'b', '-'],
        ['-', 'b', 'a', '.', '-', 'b', 'a'],
        -2.0)
    
    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.dialign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef float gapA,gapB,match,sim

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # modify matrix and traceback
    traceback[0][0] = 1
    for i in range(1,M+1):
        matrix[0][i] = matrix[0][i-1] + gopA[i-1] * scale
        traceback[0][i] = 2
    for i in range(1,N+1):
        matrix[i][0] = matrix[i-1][0] + gopB[i-1] * scale
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if proB[i-1] in r and proA[j-1] not in r and j != M:
                gapA = matrix[i-1][j] - 1000000
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]

            # calculate costs for gapB
            if proA[j-1] in r and proB[i-1] not in r and i != N:
                gapB = matrix[i][j-1] - 1000000
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]

            # calculate costs for match
            # get the score
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # check for similar prostrings
            if proA[j-1] == proB[i-1]:
                match += matrix[i-1][j-1] + match * factor
            elif proA[j-1] in r and proB[i-1] not in r:
                match += matrix[i-1][j-1] - 1000000
            elif proA[j-1] not in r and proB[i-1] in r:
                match += matrix[i-1][j-1] - 1000000
            elif abs(ord(proA[j-1])-ord(proB[i-1])) >= 2:
                match += matrix[i-1][j-1] + match * factor / 2
            else:
                match += matrix[i-1][j-1]

            # determine minimal cost
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[N][M]

    # carry out the traceback
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA += ['-']
            almB += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA += [seqA[j-1]]
            almB += [seqB[i-1]]
            i -= 1
            j -= 1
        else:
            almA += [seqA[j-1]]
            almB += ['-']
            j -= 1

    # turn alignments back
    almA,almB = almA[::-1],almB[::-1]

    # return alignments
    return almA,almB,sim

def semi_globalign(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer
        ):
    """
    Carry out semi-global alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    
    Notes
    -----
    This is the function that is called to carry out semi-global alignment
    analyses (keyword "overlap")
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. Semi-global alignment means
    that the suffixes or prefixes in one of the words are not penalized. 

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score.

    Examples
    --------
    We compare globalign with semi_globalign::
        >>> from lingpy.algorithm.cython.calign import globalign, semi_globalign
        >>> globalign(list('abab'), list('baba'), 4 * [-1], 4 * [-1], 'aaaa', 'aaaa', 4, 4, 1, 0, {("a","b"):-1, ("b","a"): -1, ("a","a"): 1, ("b", "b"): 1})
        (['a', 'b', 'a', 'b', '-'], ['-', 'b', 'a', 'b', 'a'], 1.0)
        >>> semi_globalign(list('abab'), list('baba'), 4 * [-1], 4 * [-1], 'aaaa', 'aaaa', 4, 4, 1, 0, {("a","b"):-1, ("b","a"): -1, ("a","a"): 1, ("b", "b"): 1})
        (['a', 'b', 'a', 'b', '-'], ['-', 'b', 'a', 'b', 'a'], 3.0)

    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.dialign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef float gapA,gapB,match,sim

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # modify matrix and traceback
    traceback[0][0] = 1
    for i in range(1,M+1):
        traceback[0][i] = 2
    for i in range(1,N+1):
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if j == M:
                gapA = matrix[i-1][j]
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]

            # calculate costs for gapB
            if i == N:
                gapB = matrix[i][j-1]
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]

            # calculate costs for match

            # get the score
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match += matrix[i-1][j-1] + match * factor
            elif abs(ord(proA[j-1])-ord(proB[i-1])) <= 2:
                match += matrix[i-1][j-1] + match * factor / 2
            else:
                match += matrix[i-1][j-1]

            # determine minimal cost
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[N][M]

    # carry out the traceback
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA += ['-']
            almB += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA += [seqA[j-1]]
            almB += [seqB[i-1]]
            i -= 1
            j -= 1
        else:
            almA += [seqA[j-1]]
            almB += ['-']
            j -= 1

    # turn alignments back
    almA,almB = almA[::-1],almB[::-1]

    # return alignments
    return almA,almB,sim

def secondary_semi_globalign(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer,
        str r # restricted_chars
        ):
    """
    Carry out semi-global alignment of two sequences with sensitivity to secondary sequence structures.
    
    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    r : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.

    
    Notes
    -----
    This is the function that is called to carry out semi-global alignment
    analyses (keyword "overlap")
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. Semi-global alignment means
    that the suffixes or prefixes in one of the words are not penalized. 

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score.
    
    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.dialign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef float gapA,gapB,match,sim

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # modify matrix and traceback
    traceback[0][0] = 1
    for i in range(1,M+1):
        matrix[0][i] = matrix[0][i-1] + gopA[i-1] * scale
        traceback[0][i] = 2
    for i in range(1,N+1):
        matrix[i][0] = matrix[i-1][0] + gopB[i-1] * scale
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if j == M:
                gapA = matrix[i-1][j]
            elif proB[i-1] in r and proA[j-1] not in r and j != M:
                gapA = matrix[i-1][j] - 1000000
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]

            # calculate costs for gapB
            if i == N:
                gapB = matrix[i][j-1]
            elif proA[j-1] in r and proB[i-1] not in r and i != N:
                gapB = matrix[i][j-1] - 1000000
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]

            # calculate costs for match
            # get the score
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # check for similar prostrings
            if proA[j-1] == proB[i-1]:
                match += matrix[i-1][j-1] + match * factor
            elif proA[j-1] in r and proB[i-1] not in r:
                match += matrix[i-1][j-1] - 1000000
            elif proA[j-1] not in r and proB[i-1] in r:
                match += matrix[i-1][j-1] - 1000000
            elif abs(ord(proA[j-1])-ord(proB[i-1])) <= 2:
                match += matrix[i-1][j-1] + match * factor / 2
            else:
                match += matrix[i-1][j-1]

            # determine minimal cost
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[N][M]

    # carry out the traceback
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA += ['-']
            almB += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA += [seqA[j-1]]
            almB += [seqB[i-1]]
            i -= 1
            j -= 1
        else:
            almA += [seqA[j-1]]
            almB += ['-']
            j -= 1

    # turn alignments back
    almA,almB = almA[::-1],almB[::-1]

    # return alignments
    return almA,almB,sim

def localign(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer
        ):
    """
    Carry out semi-global alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.

    
    Notes
    -----
    This is the function that is called to carry out local alignment
    analyses
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. Local alignment means
    that only the best matching substring between two sequences is returned
    (compare :evobib:`Smith1981`), also called the Smith-Waterman algorithm. 

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score. The alignments
        are each a list of suffix, alignment, and prefix.
    
    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.dialign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j,k,l

    # declare floats
    cdef float gapA,gapB,match,sim

    # declare char-character
    cdef str x

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # set similarity to zero
    sim = 0.0

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]

            # calculate costs for gapB
            if traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]

            # calculate costs for match

            # get the score
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match += matrix[i-1][j-1] + match * factor
            elif abs(ord(proA[j-1])-ord(proB[i-1])) <= 2:
                match += matrix[i-1][j-1] + match * factor / 2
            else:
                match += matrix[i-1][j-1]

            # determine minimal cost
            if gapA >= match and gapA >= gapB and gapA >= 0.0:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB and match >= 0.0:
                matrix[i][j] = match
                traceback[i][j] = 1
            elif gapB >= 0.0:
                matrix[i][j] = gapB
                traceback[i][j] = 2
            else:
                matrix[i][j] = 0.0
                traceback[i][j] = 0

            if matrix[i][j] >= sim:
                sim = matrix[i][j]
                k = i
                l = j

    # get the similarity
    sim = matrix[k][l]

    # reset i,j
    i,j = k,l

    # append stuff to almA and almB
    almA += [[x for x in seqA[j:]]]
    almB += [[x for x in seqB[i:]]]

    # append empty seq for alms to almA and almB
    almA += [[]]
    almB += [[]]
    
    while traceback[i][j] != 0:
        if traceback[i][j] == 3:
            
            almA[1] += ['-']
            almB[1] += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA[1] += [seqA[j-1]]
            almB[1] += [seqB[i-1]]
            i -= 1
            j -= 1
        
        elif traceback[i][j] == 2:
            almA[1] += [seqA[j-1]]
            almB[1] += ['-']
            j -= 1
        else:
            break
    
    # revert the alms
    almA[1] = almA[1][::-1]
    almB[1] = almB[1][::-1]

    # append the rest
    almA += [[x for x in seqA[0:j]]]
    almB += [[x for x in seqB[0:i]]]

    # return alignments
    return almA[::-1],almB[::-1],sim

def secondary_localign(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer,
        str r # restricted_chars
        ):
    """
    Carry out lobal alignment of two sequences with sensitivity to secondary sequence structures.
    
    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    r : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
        
    Notes
    -----
    This is the function that is called to carry out local alignment
    analyses
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. Local alignment means
    that only the best matching substring between two sequences is returned
    (compare :evobib:`Smith1981`), also called the Smith-Waterman algorithm. 

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score. The alignments
        are each a list of suffix, alignment, and prefix.
    
    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.dialign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j,k,l

    # declare floats
    cdef float gapA,gapB,match,sim

    # declare char-character
    cdef str x

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # set similarity to zero
    sim = 0.0

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if proB[i-1] in r and proA[j-1] not in r and j != M:
                gapA = matrix[i-1][j] - 1000000
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]

            # calculate costs for gapB
            if proA[j-1] in r and proB[i-1] not in r and j != N:
                gapB = matrix[i][j-1] - 1000000
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]

            # calculate costs for match

            # get the score
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match += matrix[i-1][j-1] + match * factor
            elif proA[j-1] in r and proB[i-1] not in r:
                match += matrix[i-1][j-1] - 1000000
            elif proA[j-1] not in r and proB[i-1] in r:
                match += matrix[i-1][j-1] - 1000000
            elif abs(ord(proA[j-1])-ord(proB[i-1])) <= 2:
                match += matrix[i-1][j-1] + match * factor / 2
            else:
                match += matrix[i-1][j-1]

            # determine minimal cost
            if gapA >= match and gapA >= gapB and gapA >= 0.0:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB and match >= 0.0:
                matrix[i][j] = match
                traceback[i][j] = 1
            elif gapB >= 0.0:
                matrix[i][j] = gapB
                traceback[i][j] = 2
            else:
                matrix[i][j] = 0.0
                traceback[i][j] = 0

            if matrix[i][j] >= sim:
                sim = matrix[i][j]
                k = i
                l = j

    # get the similarity
    sim = matrix[k][l]
    
    # reset i,j
    i,j = k,l

    # append stuff to almA and almB
    almA += [[x for x in seqA[j:]]]
    almB += [[x for x in seqB[i:]]]

    # append empty seq for alms to almA and almB
    almA += [[]]
    almB += [[]]

    while traceback[i][j] != 0:
        if traceback[i][j] == 3:
            
            almA[1] += ['-']
            almB[1] += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA[1] += [seqA[j-1]]
            almB[1] += [seqB[i-1]]
            i -= 1
            j -= 1
        
        elif traceback[i][j] == 2:
            almA[1] += [seqA[j-1]]
            almB[1] += ['-']
            j -= 1
        else:
            break

    # revert the alms
    almA[1] = almA[1][::-1]
    almB[1] = almB[1][::-1]

    # append the rest
    almA += [[x for x in seqA[0:j]]]
    almB += [[x for x in seqB[0:i]]]

    # return alignments
    return almA[::-1],almB[::-1],sim

def dialign(
        list seqA,
        list seqB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer
        ):
    """
    Carry out dialign alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
        
    Notes
    -----
    This is the function that is called to carry out local dialign alignment
    analyses (keyword "dialign")
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. Dialign (see
    :evobib:`Morgenstern1996`) is an alignment algorithm that does not require
    gap penalties and generally works in a rather local fashion.

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score.

    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.secondary_dialign

    """

    # declare integers
    cdef int i,j,k,l,o,p

    # declare floats
    cdef float gapA,gapB,match,sim,tmp_match

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # modify matrix and traceback
    traceback[0][0] = 1
    for i in range(1,M+1):
        traceback[0][i] = 2
    for i in range(1,N+1):
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            gapA = matrix[i-1][j]

            # calculate costs for gapB
            gapB = matrix[i][j-1]

            # calculate costs for match
            sim = 0.0
            o = 1
            for k in range(min(i,j)):
                match = matrix[i-k-1][j-k-1]
                for l in range(k,-1,-1):
                    # get temporary match
                    tmp_match = scorer[seqA[j-l-1],seqB[i-l-1]]
                    # check for common prostrings
                    if proA[j-l-1] == proB[i-l-1]:
                        tmp_match = tmp_match * ( 1 + factor )
                    elif abs(ord(proA[j-l-1]) - ord(proB[i-l-1])) <= 2:
                        tmp_match = tmp_match * ( 1 + factor / 2 )
                    match += tmp_match

                p = k+1
                if match > sim:
                    sim = match
                    o = p

            # determine minimal cost
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[N][M]

    # carry out the traceback
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA += ['-']
            almB += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA += [seqA[j-1]]
            almB += [seqB[i-1]]
            i -= 1
            j -= 1
        else:
            almA += [seqA[j-1]]
            almB += ['-']
            j -= 1

    # turn alignments back
    almA,almB = almA[::-1],almB[::-1]

    # return alignments
    return almA,almB,sim

def secondary_dialign(
        list seqA,
        list seqB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        object scorer,
        str r # restricted chars
        ):
    """
    Carry out dialign alignment of two sequences with sensitivity for secondary \
    sequence structures.

    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    M, N : int
        The lengths of seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`~lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    r : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
        
    Notes
    -----
    This is the function that is called to carry out local dialign alignment
    analyses (keyword "dialign")
    when using many of LingPy's classes for alignment analyses which is at the
    same time sensitive for secondary sequence structures (see the description
    of secondary alignment in :evobib:`List2014d` for details), like
    :py:class:`~lingpy.align.pairwise.Pairwise`,
    :py:class:`~lingpy.align.multiple.Multiple`, or
    :py:class:`~lingpy.compare.lexstat.LexStat`. Dialign (see
    :evobib:`Morgenstern1996`) is an alignment algorithm that does not require
    gap penalties and generally works in a rather local fashion.

    Returns
    -------
    alignment : tuple
        A tuple of the two alignments and the alignment score.

    See also
    --------
    ~lingpy.algorithm.cython.calign.globalign
    ~lingpy.algorithm.cython.calign.secondary_globalign
    ~lingpy.algorithm.cython.calign.semi_globalign
    ~lingpy.algorithm.cython.calign.secondary_semi_globalign
    ~lingpy.algorithm.cython.calign.localign
    ~lingpy.algorithm.cython.calign.secondary_localign
    ~lingpy.algorithm.cython.calign.dialign

    """

    # declare integers
    cdef int i,j,k,l,o,p

    # declare floats
    cdef float apA,gapB,match,sim,tmp_match

    # declare lists
    cdef list almA = []
    cdef list almB = []

    # create matrix and traceback
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # modify matrix and traceback
    traceback[0][0] = 1
    for i in range(1,M+1):
        traceback[0][i] = 2
    for i in range(1,N+1):
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if proB[i-1] in r and proA[j-1] not in r and j != M:
                gapA = matrix[i-1][j] - 1000000
            else:
                gapA = matrix[i-1][j]

            # calculate costs for gapB
            if proA[j-1] in r and proB[i-1] not in r and i != N:
                gapB = matrix[i][j-1] - 1000000
            else:
                gapB = matrix[i][j-1]

            # calculate costs for match
            sim = 0.0
            o = 1
            for k in range(min(i,j)):
                match = matrix[i-k-1][j-k-1]
                for l in range(k,-1,-1):
                    # get temporary match
                    tmp_match = scorer[seqA[j-l-1],seqB[i-l-1]]

                    # check for common prostrings
                    if proA[j-l-1] == proB[i-l-1]:
                        tmp_match += tmp_match * factor
                    elif proA[j-l-1] in r and proB[i-l-1] not in r:
                        tmp_match += -1000000
                    elif proA[j-l-1] not in r and proB[i-l-1] in r:
                        tmp_match += -1000000
                    elif abs(ord(proA[j-l-1]) - ord(proB[i-l-1])) <= 2:
                        tmp_match += tmp_match * factor / 2
                    
                    # get match
                    match += tmp_match

                p = k+1
                if match > sim:
                    sim = match
                    o = p

            # determine minimal cost
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[N][M]

    # carry out the traceback
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA += ['-']
            almB += [seqB[i-1]]
            i -= 1 
        elif traceback[i][j] == 1: 
            almA += [seqA[j-1]]
            almB += [seqB[i-1]]
            i -= 1
            j -= 1
        else:
            almA += [seqA[j-1]]
            almB += ['-']
            j -= 1

    # turn alignments back
    almA,almB = almA[::-1],almB[::-1]

    # return alignments
    return almA,almB,sim

def align_pair(
        list seqA,
        list seqB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int gop,
        float scale,
        float factor,
        object scorer,
        str mode,
        str restricted_chars,
        int distance = 0
        ):
    """
    Align a pair of sequences.
    
    Parameters
    ----------
    seqA, seqB : list
        The list containing the sequences.
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    mode : { "global", "local", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    restricted_chars : str
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
    distance : int (default=0)
        Select whether you want to calculate the normalized distance or the
        similarity between two strings (following :evobib:`Downey2008` for
        normalization).

    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity or distance.
    
    Notes
    -----
    This is a utility function that allows calls any of the four classical
    alignment functions (:py:class:`lingpy.algorithm.cython.calign.globalign`
    :py:class:`lingpy.algorithm.cython.calign.semi_globalign`,
    :py:class:`lingpy.algorithm.cython.calign.localign`,
    :py:class:`lingpy.algorithm.cython.calign.dialign`,) and their secondary counterparts.

    See also
    --------
    ~lingpy.algorithm.cython.calign.align_pairwise
    ~lingpy.algorithm.cython.calign.align_pairs


    """
    # define basic types
    cdef int i
    cdef list almA,almB
    cdef float sim,dist,simA,simB

    # get length of seqA,seqB
    cdef int M = len(seqA)
    cdef int N = len(seqB)

    # set up the gap costs
    gopA = [gop * gopA[i] for i in range(M)]
    gopB = [gop * gopB[i] for i in range(N)]

    # check for secondary structures
    if not set(restricted_chars).intersection(set(proA+proB)):

        # determine the mode
        if mode == "global":
            
            # carry out the alignment
            almA,almB,sim = globalign(
                    seqA,
                    seqB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer
                    )

        elif mode == "local":
            
            # carry out the alignment
            almA,almB,sim = localign(
                    seqA,
                    seqB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer
                    )

        elif mode == "overlap":
            
            # carry out the alignment
            almA,almB,sim = semi_globalign(
                    seqA,
                    seqB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer
                    )

        elif mode == "dialign":
            almA,almB,sim = dialign(
                    seqA,
                    seqB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer
                    )

    # check for secondary structures
    else:
        
        # determine the mode
        if mode == "global":
            
            # carry out the alignment
            almA,almB,sim = secondary_globalign(
                    seqA,
                    seqB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer,
                    restricted_chars
                    )

        elif mode == "local":
            
            # carry out the alignment
            almA,almB,sim = secondary_localign(
                    seqA,
                    seqB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer,
                    restricted_chars
                    )

        elif mode == "overlap":
            
            # carry out the alignment
            almA,almB,sim = secondary_semi_globalign(
                    seqA,
                    seqB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer,
                    restricted_chars
                    )

        elif mode == "dialign":
            almA,almB,sim = secondary_dialign(
                    seqA,
                    seqB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    scorer,
                    restricted_chars
                    )
    
    # calculate distance, if this is needed
    if distance > 0:
        simA = sum([(1.0 + factor) * scorer[seqA[i],seqA[i]] for i in range(M)])
        simB = sum([(1.0 + factor) * scorer[seqB[i],seqB[i]] for i in range(N)])

        dist = 1 - ( ( 2 * sim ) / ( simA + simB ) )
        if distance == 1:
            return almA,almB,dist
        else:
            return almA,almB,sim,dist
    else:
        return almA,almB,sim
    
def align_pairwise(
        list seqs,
        list gops,
        list pros,
        int gop,
        float scale,
        float factor,
        object scorer,
        str restricted_chars,
        str mode
        ):
    """
    Align a list of sequences pairwise.

    Parameters
    ----------
    seqs : list
        The list containing the sequences.
    gops : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    pros : list 
        The prosodic strings which have the same length as seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    mode : { "global", "local", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    r : str
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
    
    Returns
    -------
    alignments : list
        A list of tuples of size 4, containing the alignment, the
        similarity and the distance for each sequence pair.

    Notes
    -----
    This function computes alignments of all possible pairs passed in the list
    of sequences and is basically used in LingPy's module for multiple
    alignment analyses (:py:class:`lingpy.align.multiple`).

    See also
    --------
    ~lingpy.algorithm.cython.calign.align_pairs
    ~lingpy.algorithm.cython.calign.align_pair

    """
    # define basic stuff
    cdef list alignments = []
    cdef int lS = len(seqs)
    
    cdef int i,j,k,lenA,lenB
    cdef list almA,almB,seqA,seqB,gopA,gopB
    cdef float sim,simA,simB,dist
    cdef str proA,proB

    # get self-scores
    cdef list sims = [0.0 for i in range(lS)]
    cdef list lens = [0 for i in range(lS)]

    for i in range(lS):
        seqA = seqs[i]
        k = len(seqA)
        sim = sum([(1 + factor) * scorer[seqA[j],seqA[j]] for j in range(k)])
        lens[i] = k
        sims[i] = sim
        gops[i] = [gop * gops[i][j] for j in range(k)]
    
    # check for restricted chars in the beginning
    if not set(restricted_chars).intersection(set(''.join(pros))):

        if mode == "global":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        gopA,gopB = gops[i],gops[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = globalign(
                                seqA,
                                seqB,
                                gopA,
                                gopB,
                                proA,
                                proB,
                                lenA,
                                lenB,
                                scale,
                                factor,
                                scorer
                                )
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )

        elif mode == "local":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        gopA,gopB = gops[i],gops[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        # check for secondary structures
                        almA,almB,sim = localign(
                                seqA,
                                seqB,
                                gopA,
                                gopB,
                                proA,
                                proB,
                                lenA,
                                lenB,
                                scale,
                                factor,
                                scorer
                                ) 
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )

        elif mode == "overlap":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        gopA,gopB = gops[i],gops[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = semi_globalign(
                                seqA,
                                seqB,
                                gopA,
                                gopB,
                                proA,
                                proB,
                                lenA,
                                lenB,
                                scale,
                                factor,
                                scorer
                                )
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )

        elif mode == "dialign":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = dialign(
                               seqA,
                               seqB,
                               proA,
                               proB,
                               lenA,
                               lenB,
                               scale,
                               factor,
                               scorer
                               )
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )
    else:
        if mode == "global":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        gopA,gopB = gops[i],gops[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = secondary_globalign(
                                seqA,
                                seqB,
                                gopA,
                                gopB,
                                proA,
                                proB,
                                lenA,
                                lenB,
                                scale,
                                factor,
                                scorer,
                                restricted_chars
                                )
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )

        elif mode == "local":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        gopA,gopB = gops[i],gops[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = secondary_localign(
                                seqA,
                                seqB,
                                gopA,
                                gopB,
                                proA,
                                proB,
                                lenA,
                                lenB,
                                scale,
                                factor,
                                scorer,
                                restricted_chars
                                ) 
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )

        elif mode == "overlap":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        gopA,gopB = gops[i],gops[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = secondary_semi_globalign(
                                seqA,
                                seqB,
                                gopA,
                                gopB,
                                proA,
                                proB,
                                lenA,
                                lenB,
                                scale,
                                factor,
                                scorer,
                                restricted_chars
                                )
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )

        elif mode == "dialign":
            # start loop
            for i in range(lS):
                for j in range(lS):
                    if i < j:
                        seqA,seqB = seqs[i],seqs[j]
                        proA,proB = pros[i],pros[j]
                        simA,simB = sims[i],sims[j]
                        lenA,lenB = lens[i],lens[j]

                        almA,almB,sim = secondary_dialign(
                               seqA,
                               seqB,
                               proA,
                               proB,
                               lenA,
                               lenB,
                               scale,
                               factor,
                               scorer,
                               restricted_chars
                               )
                        
                        # get the distance
                        dist = 1 - ( 2 * sim / ( simA + simB ) )
                        
                        # append it to list
                        alignments.append(
                                (almA,almB,sim,dist)
                                )
                    elif i == j:
                        seqA = seqs[i]
                        alignments.append(
                                (seqA,seqA,sims[i],0.0)
                                )


    return alignments
    
def align_pairs(
        list seqs,
        list gops,
        list pros,
        int gop,
        float scale,
        float factor,
        object scorer,
        str mode,
        str restricted_chars,
        int distance = 0
        ):
    """
    Align multiple sequence pairs.

    Parameters
    ----------
    seqs : list
        A two-dimensional list containing one pair of sequences each. 
    gops : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    pros : list 
        The prosodic strings which have the same length as seqA and seqB.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in seqA and seqB.
    mode : { "global", "local", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    restricted_chars : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
    distance : int (default=0)
        Select whether you want to calculate the normalized distance or the
        similarity between two strings (following :evobib:`Downey2008` for
        normalization). If you set this value to 2, both distances and
        similarities will be returned.

    Returns
    -------
    alignments : list
        A list of tuples of size 3 or 4, containing the alignments, and the
        similarity or the distance (or both, if distance is set to 2).

    Notes
    -----
    This function computes alignments of all pairs passed in the list
    of sequence pairs (a two-dimensional list with two sequences each)
    and is basically used in LingPy's module for cognate detection
    (:py:class:`lingpy.compare.lexstat.LexStat`).

    See also
    --------
    ~lingpy.algorithm.cython.calign.align_pairwise
    ~lingpy.algorithm.cython.calign.align_pair
    """
    # basic defs
    cdef int i,j,M,N,lP
    cdef list seqA,seqB,almA,almB
    cdef float sim
    cdef list alignments = []

    # get basic params
    lP = len(seqs)

    # check for restricted prostrings

    # carry out alignments
    for i in range(lP):
        # get sequences
        seqA,seqB = seqs[i][0],seqs[i][1]
        
        # get length of seqs
        M,N = len(seqA),len(seqB)
        
        # get gops
        gopA = [gop * gops[i][0][j] for j in range(M)]
        gopB = [gop * gops[i][1][j] for j in range(N)]

        # get pros
        proA,proB = pros[i][0],pros[i][1]

        # check for restricted chars
        if not set(restricted_chars).intersection(proA+proB):
            if mode == "global":
                almA,almB,sim = globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer
                       )
            elif mode == "local":
                almA,almB,sim = localign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer
                       )

            elif mode == "overlap":
                almA,almB,sim = semi_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       )

            elif mode == "dialign":
                almA,almB,sim = dialign(
                       seqA,
                       seqB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       )

        else:
            if mode == "global":
                almA,almB,sim = secondary_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )
            elif mode == "local":
                almA,almB,sim = secondary_localign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

            elif mode == "overlap":
                almA,almB,sim = secondary_semi_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

            elif mode == "dialign":
                almA,almB,sim = secondary_dialign(
                       seqA,
                       seqB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

        # calculate distances if option is chose
        if distance > 0:
            simA = sum([(1.0 + factor) * scorer[seqA[i],seqA[i]] for i in range(M)])
            simB = sum([(1.0 + factor) * scorer[seqB[i],seqB[i]] for i in range(N)])

            dist = 1 - ( ( 2 * sim ) / ( simA + simB ) )
            if distance == 1:
                alignments.append((almA,almB,dist))
            else:
                alignments.append((almA,almB,sim,dist))
        else:
            alignments.append((almA,almB,sim))
    
    return alignments

# specific methods for the alignment of profiles
def align_profile(
        list profileA,
        list profileB,
        list gopA,
        list gopB,
        str proA,
        str proB,
        int gop,
        float scale,
        float factor,
        object scorer,
        str restricted_chars,
        str mode,
        float gap_weight
        ):
    """
    Align two profiles using the basic modes.

    Parameters
    ----------
    profileA, profileB : list
        Two-dimensional list for each of the profiles. 
    gopA, gopB : list
        The gap opening penalties (individual for each sequence, therefore
        passed as a list of floats or integers).
    proA, proB : str
        The prosodic strings which have the same length as profileA and profileB.
    gop : int
        The general gap opening penalty which will be used to introduce a gap
        between the two profiles.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in the two profiles.
    restricted_chars : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
        They need to be computed by computing a consensus string from all
        prosodic strings in the profile.
    mode : { "global", "local", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    gap_weight : float
        This handles the weight that is given to gaps in a column. If you set
        it to 0, for example, this means that all gaps will be ignored when
        determining the score for two columns in the profile.

    Notes
    -----
    This function computes alignments of two profiles of multiple sequences
    (see :evobib:`Durbin2002` for details on profiles)
    and is basically used in LingPy's module for multiple alignment
    (:py:class:`lingpy.align.multiple`).

    Returns
    -------
    alignment : tuple
        The aligned profiles, and the overall similarity of the profiles.

    See also
    --------
    ~lingpy.algorithm.cython.calign.score_profile
    ~lingpy.algorithm.cython.calign.swap_score_profile

    """

    # basic defs
    cdef int i,j,k,l,M,N,O,P
    cdef float sim,count
    cdef str charA,charB
    cdef list listA,listB,almA,almB
    
    M = len(profileA)
    N = len(profileB)
    O = len(profileA[0])
    P = len(profileB[0])

    cdef dict tmp_scorer = {}

    listA = [i for i in range(M)]
    listB = [i for i in range(N)]

    for i in range(M):
        for j in range(N):
            sim = 0.0
            count = 0.0
            for k in range(O):
                for l in range(P):
                    charA = profileA[i][k]
                    charB = profileB[j][l]
                    if charA != 'X' and charB != 'X':
                        sim += scorer[charA,charB]
                        count += 1.0
                    else:
                        count += gap_weight
            tmp_scorer[i,j] = sim / count

    # get the gop
    gopA = [gop * gopA[i] for i in range(M)]
    gopB = [gop * gopB[i] for i in range(N)]
    
    if not set(restricted_chars).intersection(proA+proB):
        if mode == "global":
            almA,almB,sim = globalign(
                    listA,
                    listB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    tmp_scorer
                    )
        elif mode == "overlap":
            almA,almB,sim = semi_globalign(
                    listA,
                    listB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    tmp_scorer
                    )
        elif mode == "dialign":
            almA,almB,sim = dialign(
                    listA,
                    listB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    tmp_scorer
                    )
    else:
        if mode == "global":
            almA,almB,sim = secondary_globalign(
                    listA,
                    listB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    tmp_scorer,
                    restricted_chars
                    )
        elif mode == "overlap":
            almA,almB,sim = secondary_semi_globalign(
                    listA,
                    listB,
                    gopA,
                    gopB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    tmp_scorer,
                    restricted_chars
                    )
        elif mode == "dialign":
            almA,almB,sim = secondary_dialign(
                    listA,
                    listB,
                    proA,
                    proB,
                    M,
                    N,
                    scale,
                    factor,
                    tmp_scorer,
                    restricted_chars
                    )

    return almA,almB,sim

# functions for profile scoring
def score_profile(
        list colA,
        list colB,
        object scorer,
        float gap_weight = 0.0
        ):
    """
    Basic function for the scoring of profiles.

    Parameters
    ----------
    colA, colB : list
        The two columns of a profile.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in the two profiles.
    gap_weight : float (default=0.0)
        This handles the weight that is given to gaps in a column. If you set
        it to 0, for example, this means that all gaps will be ignored when
        determining the score for two columns in the profile.
    
    Notes
    -----
    This function handles how profiles are scored.
    
    Returns
    -------
    score : float
        The score for the profile

    See also
    --------
    ~lingpy.algorithm.cython.calign.align_profile
    ~lingpy.algorithm.cython.calign.swap_score_profile

    """
    # basic definitions
    cdef int i,j
    cdef str charA,charB

    # define the initial score
    cdef float score = 0.0

    # set a counter
    cdef float counter = 0

    # iterate over all chars
    for i,charA in enumerate(colA):
        for j,charB in enumerate(colB):
            if charA != 'X' and charB != 'X':
                score += scorer[charA,charB]
                counter += 1.0
            else:
                counter += gap_weight

    return score / counter

def swap_score_profile(
        list colA,
        list colB,
        object scorer,
        float gap_weight = 0.0,
        int swap_penalty = -5
        ):
    """
    Basic function for the scoring of profiles which contain swapped sequences.

    Parameters
    ----------
    colA, colB : list
        The two columns of a profile.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in the two profiles.
    gap_weight : float (default=0.0)
        This handles the weight that is given to gaps in a column. If you set
        it to 0, for example, this means that all gaps will be ignored when
        determining the score for two columns in the profile.
    swap_penalty : int (default=-5)
        The swap penalty applied to swapped columns.
    
    Notes
    -----
    This function handles how profiles with swapped segments are scored.

    Returns
    -------
    score : float
        The score for the profile.

    See also
    --------
    ~lingpy.algorithm.cython.calign.align_profile
    ~lingpy.algorithm.cython.calign.score_profile

    """
    # basic definitions
    cdef int i,j
    cdef str charA,charB

    # define the initial score
    cdef float score = 0.0

    # set a counter
    cdef float counter = 0

    # iterate over all chars
    for i,charA in enumerate(colA):
        for j,charB in enumerate(colB):
            if charA != 'X' and charB != 'X' and charA != '+' and charB != '+':
                score += scorer[charA,charB]
                counter += 1.0
            elif charA == '+' or charB == '+':
                if charA == '+' and charB == '+':
                    score += 0.0
                    counter += 1.0
                elif charA == 'X' or charB == 'X':
                    score += swap_penalty # this is the swap cost
                    counter += 1.0
                else:
                    score += -1000000
                    counter += 1.0
            else:
                counter += gap_weight

    return score / counter

def corrdist(
        float threshold,
        list seqs,
        list gops,
        list pros,
        int gop,
        float scale,
        float factor,
        object scorer,
        str mode,
        str restricted_chars
        ):
    """
    Create a correspondence distribution for a given language pair.

    Parameters
    ----------
    threshold : float
        The threshold of sequence distance which determines whether a sequence
        pair is included or excluded from the calculation of the distribution.
    seqs : list
        The sequences passed as a two-dimensional list of sequence pairs.
    gops : list
        The gap opening penalties, passed as individual lists of penalties for each
        sequence.
    pros : list
        The list of prosodic strings for each sequence.
    gop : int
        The general gap opening penalty which will be used to introduce a gap
        between the two profiles.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    factor : float
        The factor by which matches are increased when two segments occur in
        the same prosodic position of an alignment.
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in the two profiles.
    mode : { "global", "local", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    restricted_chars : { str }
        The string containing restricted characters. Restricted characters
        occur, as a rule, in the prosodic strings, not in the normal sequence.
        They need to be computed by computing a consensus string from all
        prosodic strings in the profile.
    
    Notes
    -----
    This function is the core of the
    :py:class:`~lingpy.compare.lexstat.LexStat` function to compute
    distributions of aligned segment pairs.
    
    Returns
    -------
    results : tuple
        A dictionary containing the distribution, and the number of included
        sequences. 
    """

    # basic defs
    cdef int i,j,M,N,lP,l
    cdef list seqA,seqB,almA,almB
    cdef float sim
    cdef dict corrs = {}

    # return number of sequences considered for initial distribution
    cdef int included = 0

    # get basic params
    lP = len(seqs)

    # check for restricted prostrings

    # carry out alignments
    for i in range(lP):
        # get sequences
        seqA,seqB = seqs[i][0],seqs[i][1]
        
        # get length of seqs
        M,N = len(seqA),len(seqB)
        
        # get gops
        gopA = [gop * gops[i][0][j] for j in range(M)]
        gopB = [gop * gops[i][1][j] for j in range(N)]

        # get pros
        proA,proB = pros[i][0],pros[i][1]

        # check for restricted chars
        if not set(restricted_chars).intersection(proA+proB):
            if mode == "global":
                almA,almB,sim = globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer
                       )
            elif mode == "local":
                almA,almB,sim = localign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer
                       )
                almA = almA[1]
                almB = almB[1]


            elif mode == "overlap":
                almA,almB,sim = semi_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       )

            elif mode == "dialign":
                almA,almB,sim = dialign(
                       seqA,
                       seqB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       )

        else:
            if mode == "global":
                almA,almB,sim = secondary_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )
            elif mode == "local":
                almA,almB,sim = secondary_localign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )
                almA = almA[1]
                almB = almB[1]

            elif mode == "overlap":
                almA,almB,sim = secondary_semi_globalign(
                       seqA,
                       seqB,
                       gopA,
                       gopB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

            elif mode == "dialign":
                almA,almB,sim = secondary_dialign(
                       seqA,
                       seqB,
                       proA,
                       proB,
                       M,
                       N,
                       scale,
                       factor,
                       scorer,
                       restricted_chars
                       )

        # calculate distances
        simA = sum([(1.0 + factor) * scorer[seqA[i],seqA[i]] for i in range(M)])
        simB = sum([(1.0 + factor) * scorer[seqB[i],seqB[i]] for i in range(N)])
        dist = 1 - ( ( 2 * sim ) / ( simA + simB ) )
        
        if dist <= threshold:
            included += 1
            l = len(almA)
            for j in range(l):
                try:
                    corrs[almA[j],almB[j]] += 1
                except:
                    corrs[almA[j],almB[j]] = 1

    return corrs, included


