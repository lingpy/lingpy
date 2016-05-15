from __future__ import unicode_literals
"""
This module provides various alignment functions in an optimized version.
"""

def nw_align(
        seqA,
        seqB,
        scorer,
        gap
        ):
    """
    Align two sequences using the Needleman-Wunsch algorithm.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as list.
    scorer : dict
        A dictionary containing tuples of two segments as key and numbers as
        values.
    gap : int
        The gap penalty.
    
    Returns
    -------
    alignment : tuple
        A of the two aligned sequences, and the similarity score.

    Notes
    -----
    This function is a very straightforward implementation of the
    Needleman-Wunsch algorithm (:evobib:`Needleman1970`). We recommend to use
    the function if you want to test your own scoring dictionaries and profit
    from a fast implementation (as we use Cython, the implementation is indeed
    faster than pure Python implementations, as long as you use Python 3 and
    have Cython installed). If you want to test the NW algorithm without
    specifying a scoring dictionary, we recommend to have a look at our wrapper
    function with the same name in the :py:class:`~lingpy.align.pairwise`
    module.

    """
    
    # get the lengths of the strings
    M = len(seqA)
    N = len(seqB)

    # define general and specific integers
# [autouncomment]     cdef int i,j
# [autouncomment]     cdef int sim # stores the similarity score

    # define values for the main loop
# [autouncomment]     cdef int gapA,gapB,match,penalty # for the loop
 
    # define values for the traceback
    almA = []
    almB = []

    # create matrix and traceback
    matrix = [[0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]

    # initialize matrix and traceback
    for i in range(1,M+1):
        matrix[0][i] = matrix[0][i-1] + gap
        traceback[0][i] = 2
    for i in range(1,N+1):
        matrix[i][0] = matrix[i-1][0] + gap
        traceback[i][0] = 3
    
    # start the main loop
    for i in range(1,N+1):
        for j in range(1,M+1):
            
            # get the penalty
            match = scorer[seqA[j-1],seqB[i-1]]
            
            # get the three scores
            gapA = matrix[i-1][j] + gap
            gapB = matrix[i][j-1] + gap
            match = matrix[i-1][j-1] + match

            # evaluate the scores
            if gapA >= match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    # get the similarity
    sim = matrix[i][j] 
    
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

    return (almA[::-1],almB[::-1],sim)

def edit_dist(
        seqA,
        seqB,
        normalized
        ):
    """
    Return the edit-distance between two strings.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as list.
    normalized : bool
        Indicate whether you want the normalized or the unnormalized edit
        distance to be returned.

    Note
    ----
    This function computes the edit distance between two type objects. We
    recommend to use it if you need a fast implementation. Otherwise,
    especially, if you want to pass strings, we recommend to have a look at the
    wrapper function with the same name in the
    :py:class:`~lingpy.align.pairwise` module.

    Returns
    -------
    dist : { int, }
        Either the normalized or the unnormalized edit distance.

    """
    
    M = len(seqA)
    N = len(seqB)
# [autouncomment]     cdef int gapA,gapB,match
# [autouncomment]     cdef int i,j,sim
# [autouncomment]     cdef float dist
    
    matrix = [[0 for i in range(M+1)] for j in range(N+1)]
    
    for i in range(1,M+1):
        matrix[0][i] = i
    for i in range(1,N+1):
        matrix[i][0] = i

    for i in range(1,N+1):
        for j in range(1,M+1):
            
            if seqA[j-1] == seqB[i-1]:
                match = matrix[i-1][j-1]
            else:
                match = matrix[i-1][j-1] + 1

            gapA = matrix[i-1][j] + 1
            gapB = matrix[i][j-1] + 1

            if gapA < match and gapA < gapB:
                matrix[i][j] = gapA
            elif match <= gapB:
                matrix[i][j] = match
            else:
                matrix[i][j] = gapB

    sim = matrix[N][M]
    
    if normalized:
        dist = float(sim) / max([M,N])
        return dist

    return sim

def sw_align(
        seqA,
        seqB,
        scorer,
        gap
        ):
    """
    Align two sequences using the Smith-Waterman algorithm.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as list.
    scorer : dict
        A dictionary containing tuples of two segments as key and numbers as
        values.
    gap : int
        The gap penalty.
    
    Returns
    -------
    alignment : tuple
        A of the two aligned sequences, and the similarity score.

    Notes
    -----
    This function is a very straightforward implementation of the
    Smith-Waterman algorithm (:evobib:`Smith1981`). We recommend to use
    the function if you want to test your own scoring dictionaries and profit
    from a fast implementation (as we use Cython, the implementation is indeed
    faster than pure Python implementations, as long as you use Python 3 and
    have Cython installed). If you want to test the SW algorithm without
    specifying a scoring dictionary, we recommend to have a look at our wrapper
    function with the same name in the :py:class:`~lingpy.align.pairwise`
    module.

    """
    # basic stuff
# [autouncomment]     cdef int i,j
# [autouncomment]     cdef float gapA,gapB

    # get the lengths of the strings
    lenA = len(seqA)
    lenB = len(seqB)

# [autouncomment]     cdef str s

    # define values for the main loop
    null = 0 # constant during the loop
    imax = 1 # for the loop
    jmax = 1 # for the loop
    max_score = 0.0 # for the loo

    # define values for the traceback
    igap = 0
    jgap = 0 
    almA = [s for s in seqA]
    almB = [s for s in seqB]
    gap_char = '-' # the gap character

    # create matrix and traceback
    matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    
    # start the main loop
    for i in range(1,lenB+1):
        for j in range(1,lenA+1):
            
            # get the penalty
            penalty = scorer[seqA[j-1],seqB[i-1]]
            
            # get the three scores
            gapA = matrix[i-1][j] + gap
            gapB = matrix[i][j-1] + gap
            match = matrix[i-1][j-1] + penalty

            # evaluate the scores
            if gapA >= match and gapA >= gapB and gapA >= null:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB and match >= null:
                matrix[i][j] = match
                traceback[i][j] = 1
            elif gapB >= null:
                matrix[i][j] = gapB
                traceback[i][j] = 2
            else:
                matrix[i][j] = null
                traceback[i][j] = null

            # check for maximal score
            if matrix[i][j] >= max_score:
                imax = i
                jmax = j
                max_score = matrix[i][j]

    # get the similarity
    sim = matrix[imax][jmax]

    # start the traceback
    i,j = imax,jmax
    igap,jgap = 0,0

    while traceback[i][j] != 0:
        if traceback[i][j] == 3:
            almA.insert(j,gap_char)
            i -= 1
            jgap += 1
        elif traceback[i][j] == 1:
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:
            almB.insert(i,gap_char)
            j -= 1
            igap += 1
        else:
            break

    # return the alignment as a of prefix, alignment, and suffix
    return (
            (
                almA[0:j],
                almA[j:jmax+jgap],
                almA[jmax+jgap:]
                ),
            (
                almB[0:i],
                almB[i:imax+igap],
                almB[imax+igap:]
                ),
            sim
            )

def we_align(
        seqA,
        seqB,
        scorer,
        gap
        ):
    """
    Align two sequences using the Waterman-Eggert algorithm.

    Parameters
    ----------
    seqA, seqB : list
        The input sequences passed as a list.
    scorer : dict
        A dictionary containing tuples of two segments as key and numbers as
        values.
    gap : 
        The gap penalty.

    Notes
    -----
    This function is a very straightforward implementation of the
    Waterman-Eggert algorithm (:evobib:`Waterman1987`). We recommend to use
    the function if you want to test your own scoring dictionaries and profit
    from a fast implementation (as we use Cython, the implementation is indeed
    faster than pure Python implementations, as long as you use Python 3 and
    have Cython installed). If you want to test the WE algorithm without
    specifying a scoring dictionary, we recommend to have a look at our wrapper
    function with the same name in the :py:class:`~lingpy.align.pairwise`
    module.

    Returns
    -------
    alignments : list
        A consisting of tuples. Each tuple gives the alignment of one of
        the subsequences of the input sequences. Each contains the
        aligned part of the first, the aligned part of the second sequence, and
        the score of the alignment.

    """
    # basic defs
# [autouncomment]     cdef int lenA,lenB,i,j,null,igap,jgap
# [autouncomment]     cdef float sim,gapA,gapB,match,max_score
# [autouncomment]     cdef str gap_char
# [autouncomment]     cdef list matrix,traceback,tracer,seqA_tokens,seqB_tokens,almA,almB
    
    # get the lengths of the strings
    lenA = len(seqA)
    lenB = len(seqB)

    # define values for the main loop
    null = 0 # constant during the loop

    # define values for the traceback
    igap = 0
    jgap = 0 
    gap_char = '-' # the gap character

    # create a tracer for positions in the matrix
    tracer = [0 for i in range(lenA+1)]

    # create matrix and traceback
    matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    
    # start the main loop
    for i in range(1,lenB+1):

        # add zero to the tracer
        tracer.append(0)

        for j in range(1,lenA+1):
            
            # get the penalty
            penalty = scorer[seqA[j-1],seqB[i-1]]
            
            # get the three scores
            gapA = matrix[i-1][j] + gap
            gapB = matrix[i][j-1] + gap
            match = matrix[i-1][j-1] + penalty

            # evaluate the scores
            if gapA >= match and gapA >= gapB and gapA >= null:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB and match >= null:
                matrix[i][j] = match
                traceback[i][j] = 1
            elif gapB >= null:
                matrix[i][j] = gapB
                traceback[i][j] = 2
            else:
                matrix[i][j] = null
                traceback[i][j] = null

            # assign the value to the tracer
            tracer.append(matrix[i][j])

    
    # make of alignments
    out = []

    # start the while loop
    while True:
        
        # get the maximal value
        max_score = max(tracer)

        # if max_val is zero, break
        if max_score == 0:
            break
        
        # get the index of the maximal value of the matrix
        idx = max([i for i in range(len(tracer)) if tracer[i] == max_score])

        # convert to matrix coordinates
        i,j = idx // (lenA+1),idx - (idx // (lenA+1)) * (lenA+1)

        # store in imax and jmax
        imax,jmax = i,j

        sim = matrix[i][j]
        
        # start the traceback
        igap,jgap = 0,0

        # make values for almA and almB
        almA = [s for s in seqA]
        almB = [s for s in seqB]

        while traceback[i][j] != 0:
            if traceback[i][j] == 3:
                almA.insert(j,gap_char)
                #tracer[i * (lenA+1) + j] = 0 # set tracer to zero
                i -= 1
                jgap += 1
            elif traceback[i][j] == 1:
                #tracer[i * (lenA+1) + j] = 0 # set tracer to zero
                i -= 1
                j -= 1
            elif traceback[i][j] == 2:
                almB.insert(i,gap_char)
                #tracer[i * (lenA+1) + j] = 0 # set tracer to zero
                j -= 1
                igap += 1
            else:
                break
        
        # store values
        imin,jmin = i,j

        # change values to 0 in the tracer
        for i in range(1,lenB+1):
            for j in range(1,lenA+1):
                if imin < i <= imax or jmin < j <= jmax:
                    tracer[i * (lenA+1) + j] = 0
                    traceback[i][j] = 0
        
        # retrieve the aligned parts of the sequences
        out.append((almA[jmin:jmax+jgap],almB[imin:imax+igap],sim))

    # return the alignment as a of prefix, alignment, and suffix
    return out

def structalign(
        seqA,
        seqB,
        restricted_char = ''
        ):
    """
    Carry out a structural alignment analysis using Dijkstra's algorithm.
    
    Parameters
    ----------
    seqA,seqB : str
        The input sequences.
    restricted_chars : (default = "")
        The characters which are used to separate secondary from primary
        segments in the input sequences. Currently, the use of restricted chars
        may fail to yield an alignment.

    Notes
    -----
    Structural alignment is hereby understood as an alignment of two sequences
    whose alphabets differ. The algorithm returns all alignments with minimal
    edit distance. Edit distance in this context refers to the number of edit
    operations that are needed in order to convert one sequence into the other,
    with repeated edit operations being penalized only once.
    """
    # get basic variables
# [autouncomment]     cdef int maxScore,thisScore,newScore,fullScore
# [autouncomment]     cdef list out,queue,alm
# [autouncomment]     cdef str restA,restB
# [autouncomment]     cdef tuple residues

    # get the max score
    maxScore = max(len(seqA),len(seqB))

    # set up the queue
    queue = [
            [
                [],
                0,
                seqA,
                seqB
                ]
            ]
    
    out = []

    # while loop
    while queue:
        
        # get the first element of the queue
        alm,thisScore,restA,restB = queue.pop(0)

        if not restA and not restB and thisScore <= maxScore:
            out += [(''.join([a[0] for a in alm]),''.join([a[1] for a in alm]))]

        # start adding match
        if restA and restB:
            residues = (restA[0],restB[0])
            if residues != (" "," ") and restricted_char in residues:
                pass
            else:
                if residues not in alm:
                    newScore = thisScore + 1
                else:
                    newScore = thisScore
                fullScore = newScore + max(len(restA)-1,len(restB)-1)

                # check for better score
                if fullScore < maxScore:
                    maxScore = fullScore

                # 
                if newScore <= maxScore:
                    queue += [[alm+[residues],newScore,restA[1:],restB[1:]]]

        # start adding gap
        if restA:
            residues = (restA[0],'-')
            if restA[0] == " " and restB and restB != seqB:
                pass
            else:
                if residues not in alm:
                    newScore = thisScore + 1
                else:
                    newScore = thisScore

                fullScore = newScore + max(len(restA)-1,len(restB))

                # check for better score
                if fullScore < maxScore:
                    maxScore = fullScore

                if newScore <= maxScore:
                    queue += [[alm+[residues],newScore,restA[1:],restB]]
        
        # add gap in a
        if restB:
            residues = ('-',restB[0])
            if restB[0] == " " and restA and restA != seqA:
                pass
            else:
                if residues not in alm:
                    newScore = thisScore + 1
                else:
                    newScore = thisScore

                fullScore = newScore + max(len(restA),len(restB)-1)

                # check for better score
                if fullScore < maxScore:
                    maxScore = fullScore
                if newScore <= maxScore:
                    queue += [[alm+[residues],newScore,restA,restB[1:]]]

    return out,maxScore

def restricted_edit_dist(
        seqA,
        seqB,
        resA,
        resB,
        normalized
        ):
    r"""
    Return the restricted edit-distance between two strings.
    
    Parameters
    ----------
    seqA, seqB : list
        The two sequences passed as list.
    resA, resB : str
        The restrictions passed as a string with the same length as the
        corresponding sequence. We note a restriction if the
        strings show different symbols in their restriction string. If the
        symbols are identical, it is modeled as a non-restriction.
    normalized : bool
        Determine whether you want to return the normalized or the unnormalized
        edit distance.

    Notes
    -----
    Restrictions follow the definition of :evobib:`Heeringa2006`: Segments that
    are not allowed to match are given a penalty of :math:`\infty`. We model
    restrictions as strings, for example consisting of letters "c" and "v". So
    the sequence "woldemort" could be modeled as "cvccvcvcc", and when aligning
    it with the sequence "walter" and its restriction string "cvccvc", the
    matching of those segments in the sequences in which the segments of the
    restriction string differ, would be heavily penalized, thus prohibiting an
    alignment of "vowels" and "consonants" ("v" and "c").
    """
    
    M = len(seqA)
    N = len(seqB)
# [autouncomment]     cdef int gapA,gapB,match
# [autouncomment]     cdef int i,j,sim
# [autouncomment]     cdef float dist
    
    # define alignments
    almA = []
    almB = []
    
    # create matrix and traceback
    matrix = [[0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]   
    
    for i in range(1,M+1):
        matrix[0][i] = i
        traceback[0][i] = 2
    for i in range(1,N+1):
        matrix[i][0] = i
        traceback[i][0] = 3

    for i in range(1,N+1):
        for j in range(1,M+1):
            
            if seqA[j-1] == seqB[i-1]:
                match = matrix[i-1][j-1]
            elif resA[j-1] == resB[i-1]:
                match = matrix[i-1][j-1] + 1
            else:
                match = matrix[i-1][j-1] + 1000

            gapA = matrix[i-1][j] + 1
            gapB = matrix[i][j-1] + 1

            if gapA < match and gapA < gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match <= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    sim = matrix[N][M]
    
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
    
    if normalized:
        dist = float(sim) / len(almA)
        return dist

    return sim


