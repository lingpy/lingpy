"""
Basic module for alignment analyses.
"""
__author__ = "Johann-Mattis List"
__data__="2012-11-26"

def _get_matrix(
        tuple gopA,
        tuple gopB,
        int M,
        int N,
        str mode,
        int scale
        ):
    """
    Initialize traceback and matrix.
    """
    cdef int i,j
    
    cdef list matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]
    
    if mode != 'local':
        traceback[0][0] = 1
        if mode == 'global':
            for i in range(1,M+1):
                matrix[0][i] = matrix[0][i-1] + gopA[i-1] * scale
                traceback[0][i] = 2
            for i in range(1,N+1):
                matrix[i][0] = matrix[i-1][0] + gopB[i-1] * scale
                traceback[i][0] = 3
        else:
            for i in range(1,M+1):
                traceback[0][i] = 2
            for i in range(1,N+1):
                traceback[i][0] = 3

    return (matrix,traceback)

def _get_global_alignments(
        int i,
        int j,
        list traceback,
        tuple seqA,
        tuple seqB
        ):
    """
    Return the alignments for a given traceback.
    """
    
    cdef list almA = []
    cdef list almB = []
    
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

    return (almA[::-1],almB[::-1])

def _get_local_alignments(
        int i,
        int j,
        list traceback,
        tuple seqA,
        tuple seqB
        ):
    """
    Return the local alignments for a given traceback.
    """

    cdef list almA = []
    cdef list almB = []
    cdef str x

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

    return (almA[::-1],almB[::-1])

def _score(
        str charA,
        str charB,
        dict scorer
        ):
    """
    Basic function that scores the comparison of chars.
    """
    return scorer[charA,charB]

def _score_simple(
        str charA,
        str charB,
        dict scorer
        ):
    if charA == charB:
        return 1.0
    else: 
        return -1.0

def _self_score(
        tuple seq,
        int M,
        dict scorer,
        factor
        ):
    """
    Return the self-score for an alignment of a string with itself.
    """
    
    cdef float s
    cdef int i
    cdef float score = 0.0

    for i in range(M):
        s = _score(seq[i],seq[i],scorer)
        score += s + s * factor

    return score

def align(
        tuple seqA,
        tuple seqB,
        tuple gopA,
        tuple gopB,
        tuple proA,
        tuple proB,
        int gop,
        float scale,
        float factor,
        dict scorer,
        str res, # restricted chars
        str mode,
        bint distance = False
        ):
    """
    Calculate alignment using basic sound class approach.

    Parameters
    ---------
    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]-
    seqA,seqB : tuple
        The input sequences that are represented as a tuple of sound-class
        string, prosodic string, and gap weights.
    gop : int
        Gap penalty.
    scale : float
        The scale for consecutive gaps.
    factor : float
        The factor for segments sharing the same prosodic characer.
    scorer : dict
        The scoring dictionary.
    res : str
        The string storing the restricted characters.
    mode : str {"global","local","overlap","dialign"}
        The string storing the mode.

    Returns
    -------
    almA,almB,sim : list,list,float
        Two lists containing the alignment, and the alignment score.
    """
    
    # declare integer variables
    cdef int i,j,k,l,o,p
    # declare floats
    cdef float gapA,gapB,match,sim
    # declare lists
    cdef list almA,almB
    
    # declare arrays
    cdef list matrix,traceback 

    # check for correct mode
    cdef tuple modes = tuple(['global','local','dialign','overlap'])
    if mode not in modes:
        raise ValueError(
                '[!] Alignment mode "{0}" is not available!'.format(
                    mode
                    )
                )

    # get the length of the sequences
    cdef int M = len(seqA)
    cdef int N = len(seqB)
    
    # get the gops
    if gopA and gopB:
        gopA = tuple([gop * gopA[i] for i in range(M)])
        gopB = tuple([gop * gopB[i] for i in range(N)])
    elif not gopA and not gopB:
        gopA = tuple(M * [gop])
        gopB = tuple(N * [gop])
    else:
        raise ValueError(
                '[!] Gap weights have been specified incorrectly!'
                )
    
    # get the pros
    if not proA:
        proA = seqA
    if not proB:
        proB = seqB

    # get matrix and traceback
    matrix,traceback = _get_matrix(
        gopA,
        gopB,
        M,
        N,
        mode,
        scale
        )

    # set sim to 0 if mode is local
    if mode == 'local':
        sim = 0.0
    
    if scorer:
        score = _score
    else:
        score = _score_simple

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):
            
            # calculate costs for gapA
            if j == M and mode == 'overlap':
                gapA = matrix[i-1][j]
            elif proB[i-1] in res and proA[j-1] not in res and j != M:
                gapA = matrix[i-1][j] - 1000000000
            elif mode == 'dialign':
                gapA = matrix[i-1][j]
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]
            
            # calculate costs for gapB
            if i == N and mode == 'overlap':
                gapB = matrix[i][j-1]
            elif proA[j-1] in res and proB[i-1] not in res and i != N:
                gapB = matrix[i][j-1] - 1000000000
            elif mode == 'dialign':
                gapB = matrix[i][j-1]
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]
           
            # calculate costs for match
            if mode != 'dialign':
                match = score(seqA[j-1],seqB[i-1],scorer)
            else:
                sim = 0.0
                o = 1
                for k in range(min(i,j)):
                    match = matrix[i-k-1][j-k-1]
                    for l in range(k,-1,-1):
                        match += score(seqA[j-1],seqB[i-1],scorer)
                    p = k+1
                    if match > sim:
                        sim = match
                        o = p

            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match = matrix[i-1][j-1] + match + match * factor
            else:
                match = matrix[i-1][j-1] + match

            # global modes
            if mode != 'local':
                if gapA > match and gapA >= gapB:
                    matrix[i][j] = gapA
                    traceback[i][j] = 3
                elif match >= gapB:
                    matrix[i][j] = match
                    traceback[i][j] = 1
                else:
                    matrix[i][j] = gapB
                    traceback[i][j] = 2
            else:
                if gapA >= match and gapA >= gapB and gapA >= 0:
                    matrix[i][j] = gapA
                    traceback[i][j] = 3
                elif match >= gapB and match >= 0:
                    matrix[i][j] = match
                    traceback[i][j] = 1
                elif gapB >= 0:
                    matrix[i][j] = gapB
                    traceback[i][j] = 2
                else:
                    matrix[i][j] = 0.0
                    traceback[i][j] = 0

                if matrix[i][j] >= sim:
                    sim = matrix[i][j]
                    k = i
                    l = j
    
    # carry out the traceback
    if mode != "local":
        sim = matrix[N][M]
        almA,almB = _get_global_alignments(
                i,
                j,
                traceback,
                seqA,
                seqB,
                )
    else:
        sim = matrix[k][l]
        almA,almB = _get_local_alignments(
                k,
                l,
                traceback,
                seqA,
                seqB
                )
    
    if not distance:
        return (almA,almB,sim)
    else:
        # calculate distance using Downey's formula
        gapA = _self_score(
                seqA,
                M,
                scorer,
                factor
                )
        gapB = _self_score(
                seqB,
                N,
                scorer,
                factor
                )

        return (almA,almB,1 - (2 * sim / (gapA+gapB)))

def nw_align(
        tuple seqA,
        tuple seqB,
        dict scorer,
        int gap = -1
        ):
    """
    Align two sequences using the Needleman-Wunsch algorithm.
    """
    
    # get the lengths of the strings
    cdef int M = len(seqA)
    cdef int N = len(seqB)

    # define lists for tokens (in case no scoring function is provided)
    #cdef list seqA_tokens,seqB_tokens
    #cdef str tA,tB

    # define general and specific integers
    cdef int i,j
    cdef int sim # stores the similarity score

    # define values for the main loop
    cdef int gapA,gapB,match,penalty # for the loop
 
    # define values for the traceback
    cdef list almA 
    cdef list almB
    #cdef str gap_char = '-' # the gap character
    
    # define lists for tokens (in case no scoring function is provided)
    cdef list seqA_tokens,seqB_tokens
    cdef str tA,tB
    
    # create scorer
    cdef dict this_scorer = {}
    seqA_tokens = list(set(seqA))
    seqB_tokens = list(set(seqB))
    if not scorer:
        for tA in seqA_tokens:
            for tB in seqB_tokens:
                if tA == tB:
                    this_scorer[tA,tB] = 1
                else:
                    this_scorer[tA,tB] = -1
    else:
        for tA in seqA_tokens:
            for tB in seqB_tokens:
                this_scorer[tA,tB] = scorer[tA,tB]

    # create matrix and traceback
    cdef list matrix = [[0 for i in range(M+1)] for j in range(N+1)]
    cdef list traceback = [[0 for i in range(M+1)] for j in range(N+1)]

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
            match = this_scorer[seqA[j-1],seqB[i-1]]
            
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

    # get the traceback
    almA,almB = _get_global_alignments(
            i,
            j,
            traceback,
            seqA,
            seqB
            )

    # return the alignment as a tuple of prefix, alignment, and suffix
    return (almA,almB,sim)


#def sw_align(
#        tuple seqA,
#        tuple seqB,
#        dict scorer,
#        int gop = -1,
#        int gep = -1
#        ):
#    """
#    Align two sequences using the Smith-Waterman algorithm.
#    """
#    
#    # get the lengths of the strings
#    cdef int lenA = len(seqA)
#    cdef int lenB = len(seqB)
#
#    # define lists for tokens (in case no scoring function is provided)
#    cdef list seqA_tokens,seqB_tokens
#    cdef str tA,tB
#
#    # define general and specific integers
#    cdef int i,j
#    cdef int sim # stores the similarity score
#
#    # define values for the main loop
#    cdef int gapA,gapB,match # for the loop
#    cdef int null = 0 # constant during the loop
#    cdef int imax = 1 # for the loop
#    cdef int jmax = 1 # for the loop
#    cdef int max_score = 0 # for the loo
#
#    # define values for the traceback
#    cdef int igap = 0
#    cdef int jgap = 0 
#    cdef list almA = seqA[:]
#    cdef list almB = seqB[:] 
#    cdef str gap_char = '-' # the gap character
#
#    # create matrix and traceback
#    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
#    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
#
#    # create scorer, if it is empty
#    if not scorer:
#        seqA_tokens = list(set(seqA))
#        seqB_tokens = list(set(seqB))
#        for tA in seqA_tokens:
#            for tB in seqB_tokens:
#                if tA == tB:
#                    scorer[tA,tB] = 1
#                else:
#                    scorer[tA,tB] = -1
#    
#    # start the main loop
#    for i in range(1,lenB+1):
#        for j in range(1,lenA+1):
#            
#            # get the penalty
#            match = scorer[seqA[j-1],seqB[i-1]]
#            
#            # get the three scores
#            gapA = matrix[i-1][j] + gap
#            gapB = matrix[i][j-1] + gap
#            match = matrix[i-1][j-1] + penalty
#
#            # evaluate the scores
#            if gapA >= match and gapA >= gapB and gapA >= null:
#                matrix[i][j] = gapA
#                traceback[i][j] = 3
#            elif match >= gapB and match >= null:
#                matrix[i][j] = match
#                traceback[i][j] = 1
#            elif gapB >= null:
#                matrix[i][j] = gapB
#                traceback[i][j] = 2
#            else:
#                matrix[i][j] = null
#                traceback[i][j] = null
#
#            # check for maximal score
#            if matrix[i][j] >= max_score:
#                imax = i
#                jmax = j
#                max_score = matrix[i][j]
#
#    # get the similarity
#    sim = matrix[imax][jmax]
#
#    # start the traceback
#    i,j = imax,jmax
#    igap,jgap = 0,0
#
#    while traceback[i][j] != 0:
#        if traceback[i][j] == 3:
#            almA.insert(j,gap_char)
#            i -= 1
#            jgap += 1
#        elif traceback[i][j] == 1:
#            i -= 1
#            j -= 1
#        elif traceback[i][j] == 2:
#            almB.insert(i,gap_char)
#            j -= 1
#            igap += 1
#        else:
#            break
#
#    # return the alignment as a tuple of prefix, alignment, and suffix
#    return (
#            (
#                almA[0:j],
#                almA[j:jmax+jgap],
#                almA[jmax+jgap:]
#                ),
#            (
#                almB[0:i],
#                almB[i:imax+igap],
#                almB[imax+igap:]
#                ),
#            sim
#            )

def edit_dist(
        tuple seqA,
        tuple seqB,
        bint normalized = False
        ):
    """
    Return the edit-distance between two strings.
    """
    
    cdef int M = len(seqA)
    cdef int N = len(seqB)
    cdef int gapA,gapB,match
    cdef int i,j,sim
    cdef float dist
    
    cdef list matrix = [[0 for i in range(M+1)] for j in range(N+1)]
    
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
        dist = sim / float(max([M,N]))
        return dist

    return sim

