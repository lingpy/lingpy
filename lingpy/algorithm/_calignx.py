"""
Basic module for alignment analyses.
"""
__author__ = "Johann-Mattis List"
__data__="2012-11-26"

def _get_matrix(
        gopA,
        gopB,
        M,
        N,
        mode,
        scale
        ):
    """
    Initialize traceback and matrix.
    """
#     i,j
    
    matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]
    
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
        i,
        j,
        traceback,
        seqA,
        seqB
        ):
    """
    Return the alignments for a given traceback.
    """
    
    almA = []
    almB = []
    
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
        i,
        j,
        traceback,
        seqA,
        seqB
        ):
    """
    Return the local alignments for a given traceback.
    """

    almA = []
    almB = []
#     x

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
        charA,
        charB,
        scorer
        ):
    """
    Basic function that scores the comparison of chars.
    """
    return scorer[charA,charB]

def _score_simple(
        charA,
        charB,
        scorer
        ):
    if charA == charB:
        return 1.0
    else: 
        return -1.0

def _score_profile(
        colA,
        colB,
        scorer,
        gap_weight = 0.0
        ):
    """
    Basic function for the scoring of profiles.
    """
    # basic definitions
#     i,j
#     charA,charB

    # define the initial score
    score = 0.0

    # set a counter
    counter = 0

    # iterate over all chars
    for i,charA in enumerate(colA):
        for j,charB in enumerate(colB):
            if charA != 'X' and charB != 'X':
                score += scorer[charA,charB]
                counter += 1.0
            else:
                counter += gap_weight

    return score / counter

def _swap_score_profile(
        colA,
        colB,
        scorer,
        gap_weight = 0.0,
        swap_penalty = -5
        ):
    """
    Basic function for the scoring of profiles.
    """
    # basic definitions
#     i,j
#     charA,charB

    # define the initial score
    score = 0.0

    # set a counter
    counter = 0

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

def _basic_score_profile(
        colA,
        colB,
        scorer,
        gop,
        gap_weight = 0.0
        ):
    """
    Basic function for the scoring of profiles.
    """
    # basic definitions
#     i,j
#     charA,charB

    # define the initial score
    score = 0.0

    # set a counter
    counter = 0

    # iterate over all chars
    for i,charA in enumerate(colA):
        for j,charB in enumerate(colB):
            if charA != 'X' and charB != 'X':
                score += scorer[charA,charB]
                counter += 1.0
            elif charA == 'X' and charB == 'X':
                counter += gap_weight
            else:
                score += float(gop)
                counter += 1.0
    return score / counter

def _self_sc_score(
        seq,
        M,
        scorer,
        factor
        ):
    """
    Return the self-score for an alignment of a string with itself.
    """
    
#     s
#     i
    score = 0.0

    for i in range(M):
        s = _score(seq[i],seq[i],scorer)
        score += s + s * factor

    return score

def _self_basic_score(
        seq,
        M,
        scorer,
        ):
    """
    Return the self-score for an basic alignment with itself.
    """
    
#     s
#     i
    score = 0.0

    # check for correct scorer
    if not scorer:
        scoref = _score_simple
    else:
        scoref = _score

    for i in range(M):
        s = scoref(seq[i],seq[i],scorer)
        score += s

    return score


def sc_align(
        seqA, # sequence
        seqB, # sequence
        gopA, # gap score modifier
        gopB, # gap score modifier
        proA, # prosodic string
        proB, # prosodic string
        gop,
        scale,
        factor,
        scorer,
        res, # restricted chars
        mode,
        distance = False,
        both = False # return both similarity and distance
        ):
    """
    Calculate alignment using basic sound class approach.

    Parameters
    ---------
    seqA,seqB : tuple
        The input sequences that are represented as a of sound-class
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
    mode : {"global","local","overlap","dialign"}
        The string storing the mode.

    Returns
    -------
    almA,almB,sim : list,list,float
        Two lists containing the alignment, and the alignment score.
    """
    
    # declare integer variables
#     i,j,k,l,o,p
    # declare floats
#     gapA,gapB,match,sim
    # declare lists
#     almA,almB
    
    # declare arrays
#     matrix,traceback 

    # check for correct mode
    modes = ['global','local','dialign','overlap']
    if mode not in modes:
        raise ValueError(
                '[!] Alignment mode "{0}" is not available!'.format(
                    mode
                    )
                )

    # get the length of the sequences
    M = len(seqA)
    N = len(seqB)
    
    # get the gops
    gopA = [gop * gopA[i] for i in range(M)]
    gopB = [gop * gopB[i] for i in range(N)]

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
                match = _score(seqA[j-1],seqB[i-1],scorer)
            else:
                sim = 0.0
                o = 1
                for k in range(min(i,j)):
                    match = matrix[i-k-1][j-k-1]
                    for l in range(k,-1,-1):
                        match += _score(seqA[j-1],seqB[i-1],scorer)
                    p = k+1
                    if match > sim:
                        sim = match
                        o = p

            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match = matrix[i-1][j-1] + match + match * factor
            elif abs(ord(proA[j-1])-ord(proB[i-1])) >= 2:
                match = matrix[i-1][j-1] + match + match * factor / 2
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
        gapA = _self_sc_score(
                seqA,
                M,
                scorer,
                factor
                )
        gapB = _self_sc_score(
                seqB,
                N,
                scorer,
                factor
                )
        if both:
            return (almA,almB,sim, 1 - (2 * sim / (gapA+gapB)))
        return (almA,almB,1 - (2 * sim / (gapA+gapB)))

def align_sequences_pairwise(
        seqs,
        gops,
        pros,
        gop,
        scale,
        factor,
        scorer,
        res,
        mode
        ):
    """
    Align sequences pairwise and return their alignment score.
    """
    alignments = []
    lS = len(seqs)
    
#     i,j
#     almA,almB,seqA,seqB,gopA,gopB
#     sim,simA,simB,dist
#     proA,proB

    # get the self-scores of the sequences
    self_scores = []
    for i in range(lS):
        seqA = seqs[i]
        sim = _self_sc_score(
                seqA,
                len(seqA),
                scorer,
                factor
                )
        self_scores.append(sim)

    for i in range(lS):
        for j in range(lS):
            if i < j:
                seqA,seqB = seqs[i],seqs[j]
                gopA,gopB = gops[i],gops[j]
                proA,proB = pros[i],pros[j]

                # get the self-scores
                simA,simB = self_scores[i],self_scores[j]

                # carry out alignment
                almA,almB,sim = sc_align(
                        seqA,
                        seqB,
                        gopA,
                        gopB,
                        proA,
                        proB,
                        gop,
                        scale,
                        factor,
                        scorer,
                        res,
                        mode
                        )

                # get the distance
                dist = 1 - (2 * sim / (simA + simB))

                # append everything to list
                alignments.append(
                        (almA,almB,sim,dist)
                        )
            elif i == j:
                seqA = seqs[i]
                alignments.append(
                        (seqA,seqA,self_scores[i],0.0)
                        )

    return alignments

def profile_align(
        seqA, # sequence
        seqB, # sequence
        gopA, # gap score modifier
        gopB, # gap score modifier
        proA, # prosodic string
        proB, # prosodic string
        gop,
        scale,
        factor,
        scorer,
        res, # restricted chars
        mode,
        gap_weight
        ):
    """
    Calculate alignment for profiles.

    Parameters
    ---------
    seqA,seqB : tuple
        The input sequences that are represented as a of sound-class
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
    mode : {"global","local","overlap","dialign"}
        The string storing the mode.

    Returns
    -------
    almA,almB,sim : list,list,float
        Two lists containing the alignment, and the alignment score.
    """
    
    # declare integer variables
#     i,j,k,l,o,p
    # declare floats
#     gapA,gapB,match,sim
    # declare lists
#     almA,almB
    
    # declare arrays
#     matrix,traceback 

    # check for correct mode
    modes = tuple(['global','dialign','overlap'])
    if mode not in modes:
        raise ValueError(
                '[!] Alignment mode "{0}" is not available!'.format(
                    mode
                    )
                )

    # get the length of the sequences
    M = len(seqA)
    N = len(seqB)
    
    # get the gops
    gopA = [gop * gopA[i] for i in range(M)]
    gopB = [gop * gopB[i] for i in range(N)]

    # get matrix and traceback
    matrix,traceback = _get_matrix(
        gopA,
        gopB,
        M,
        N,
        mode,
        scale
        )

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
                match = _score_profile(
                        seqA[j-1],
                        seqB[i-1],
                        scorer,
                        gap_weight
                        )
            else:
                sim = 0.0
                o = 1
                for k in range(min(i,j)):
                    match = matrix[i-k-1][j-k-1]
                    for l in range(k,-1,-1):
                        match += _score_profile(
                                seqA[j-1],
                                seqB[i-1],
                                scorer,
                                gap_weight
                                )
                    p = k+1
                    if match > sim:
                        sim = match
                        o = p

            # check for similar prostring
            if proA[j-1] == proB[i-1]:
                match = matrix[i-1][j-1] + match + match * factor
            elif abs(ord(proA[j-1])-ord(proB[i-1])) >= 2:
                match = matrix[i-1][j-1] + match + match * factor * 0.5
            else:
                match = matrix[i-1][j-1] + match

            # global modes
            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2
    
    # carry out the traceback
    sim = matrix[N][M]
    almA,almB = _get_global_alignments(
            i,
            j,
            traceback,
            seqA,
            seqB,
            )
    
    return (almA,almB,sim)

def basic_profile_align(
        seqA, # sequence
        seqB, # sequence
        gop,
        scale,
        scorer,
        mode,
        gap_weight,
        gopA = [],
        gopB = []
        ):
    """
    Calculate alignment for profiles.

    Parameters
    ---------
    seqA,seqB : tuple
        The input sequences that are represented as a of sound-class
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
    mode : {"global","local","overlap","dialign"}
        The string storing the mode.

    Returns
    -------
    almA,almB,sim : list,list,float
        Two lists containing the alignment, and the alignment score.
    """
    
    # declare integer variables
#     i,j,k,l,o,p
    # declare floats
#     gapA,gapB,match,sim
    # declare lists
#     almA,almB
    
    # declare arrays
#     matrix,traceback 

    # check for correct mode
    modes = ['global','dialign','overlap']
    if mode not in modes:
        raise ValueError(
                '[!] Alignment mode "{0}" is not available!'.format(
                    mode
                    )
                )

    # get the length of the sequences
    M = len(seqA)
    N = len(seqB)
    
    # get the gops
    if not gopA and not gopB:
        gopA = [gop for i in range(M)]
        gopB = [gop for i in range(N)]

    # get matrix and traceback
    matrix,traceback = _get_matrix(
        gopA,
        gopB,
        M,
        N,
        mode,
        scale
        )

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):
            
            # calculate costs for gapA
            if j == M and mode == 'overlap':
                gapA = matrix[i-1][j]
            elif mode == 'dialign':
                gapA = matrix[i-1][j]
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + _basic_score_profile(
                        ['X'],
                        seqA[j-1],
                        scorer,
                        gopB[i-1],
                        gap_weight
                        ) * scale
            else:
                gapA = matrix[i-1][j] + _basic_score_profile(
                        ['X'],
                        seqA[j-1],
                        scorer,
                        gopB[i-1],
                        gap_weight
                        )
            
            # calculate costs for gapB
            if i == N and mode == 'overlap':
                gapB = matrix[i][j-1]
            elif mode == 'dialign':
                gapB = matrix[i][j-1]
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + _basic_score_profile(
                        ['X'],
                        seqB[i-1],
                        scorer,
                        gopA[j-1],
                        gap_weight ) * scale
            else:
                gapB = matrix[i][j-1] + _basic_score_profile(
                        ['X'],
                        seqB[i-1],
                        scorer,
                        gopA[j-1],
                        gap_weight
                        )
           
            # calculate costs for match
            if mode != 'dialign':
                match = matrix[i-1][j-1]
                match += _basic_score_profile(
                        seqA[j-1],
                        seqB[i-1],
                        scorer,
                        gop,
                        gap_weight
                        )
            else:
                sim = 0.0
                o = 1
                for k in range(min(i,j)):
                    match = matrix[i-k-1][j-k-1]
                    for l in range(k,-1,-1):
                        match += _basic_score_profile(
                                seqA[j-1],
                                seqB[i-1],
                                scorer,
                                gop,
                                gap_weight
                                )
                    p = k+1
                    if match > sim:
                        sim = match
                        o = p

            # global modes
            if gapA > match and gapA > gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    
    # carry out the traceback
    sim = matrix[N][M]
    almA,almB = _get_global_alignments(
            i,
            j,
            traceback,
            seqA,
            seqB,
            )
    
    return (almA,almB,sim)


def basic_align(
        seqA, # sequence
        seqB, # sequence
        gop,
        scale,
        scorer,
        mode,
        distance = False,
        both = False,
        gopA = [],
        gopB = []
        ):
    """
    Calculate alignment using simple character approach.

    Parameters
    ---------
    seqA,seqB : tuple
        The input sequences that are represented as a of sound-class
        string, prosodic string, and gap weights.
    gop : int
        Gap penalty.
    scale : float
        The scale for consecutive gaps.
    scorer : dict
        The scoring dictionary.
    res : str
        The string storing the restricted characters.
    mode : {"global","local","overlap","dialign"}
        The string storing the mode.

    Returns
    -------
    almA,almB,sim : list,list,float
        Two lists containing the alignment, and the alignment score.
    """
    
    # declare integer variables
#     i,j,k,l,o,p
    # declare floats
#     gapA,gapB,match,sim
    # declare lists
#     almA,almB
    
    # declare arrays
#     matrix,traceback 

    # check for correct mode
    modes = ['global','local','dialign','overlap']
    if mode not in modes:
        raise ValueError(
                '[!] Alignment mode "{0}" is not available!'.format(
                    mode
                    )
                )

    # get the length of the sequences
    M = len(seqA)
    N = len(seqB)
    
    if not gopA and not gopB:
        gopA = M * [gop]
        gopB = N * [gop]

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
            elif mode == 'dialign':
                gapA = matrix[i-1][j]
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gopB[i-1] * scale
            else:
                gapA = matrix[i-1][j] + gopB[i-1]
            
            # calculate costs for gapB
            if i == N and mode == 'overlap':
                gapB = matrix[i][j-1]
            elif mode == 'dialign':
                gapB = matrix[i][j-1]
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gopA[j-1] * scale
            else:
                gapB = matrix[i][j-1] + gopA[j-1]
           
            # calculate costs for match
            if mode != 'dialign':
                match = matrix[i-1][j-1] + score(seqA[j-1],seqB[i-1],scorer)
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
        gapA = _self_basic_score(
                seqA,
                M,
                scorer
                )
        gapB = _self_basic_score(
                seqB,
                N,
                scorer
                )

        if both:
            return (almA,almB,sim,1 - (2 * sim / (gapA+gapB)))
        else:
            return (almA,almB,1 - (2 * sim / (gapA+gapB)))

def nw_align(
        seqA,
        seqB,
        scorer,
        gap = -1
        ):
    """
    Align two sequences using the Needleman-Wunsch algorithm.
    """
    
    # get the lengths of the strings
    M = len(seqA)
    N = len(seqB)

    # define lists for tokens (in case no scoring function is provided)
#     #seqA_tokens,seqB_tokens
#     #tA,tB

    # define general and specific integers
#     i,j
#     sim # stores the similarity score

    # define values for the main loop
#     gapA,gapB,match,penalty # for the loop
 
    # define values for the traceback
#     almA 
#     almB
    #gap_char = '-' # the gap character
    
    # define lists for tokens (in case no scoring function is provided)
#     seqA_tokens,seqB_tokens
#     tA,tB
    
    # create scorer
    this_scorer = {}
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

    # return the alignment as a of prefix, alignment, and suffix
    return (almA,almB,sim)

def edit_dist(
        seqA,
        seqB,
        normalized = False
        ):
    """
    Return the edit-distance between two strings.
    """
    
    M = len(seqA)
    N = len(seqB)
#     gapA,gapB,match
#     i,j,sim
#     dist
    
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
        dist = sim / max([M,N])
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
    """
    # basic stuff
#     i,j
#     gapA,gapB

    # get the lengths of the strings
    lenA = len(seqA)
    lenB = len(seqB)

    # define lists for tokens (in case no scoring function is provided)

    # define general and specific integers

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

    # create scorer, if it is empty
    if not scorer:
        seqA_tokens = list(set(seqA))
        seqB_tokens = list(set(seqB))
        for tA in seqA_tokens:
            for tB in seqB_tokens:
                if tA == tB:
                    scorer[tA,tB] = 1
                else:
                    scorer[tA,tB] = -1
    
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
    """
    # basic defs
#     lenA,lenB,i,j,null,igap,jgap
#     sim,gapA,gapB,match,max_score
#     gap_char
#     matrix,traceback,tracer,seqA_tokens,seqB_tokens,almA,almB
    
    # get the lengths of the strings
    lenA = len(seqA)
    lenB = len(seqB)

    # define lists for tokens (in case no scoring function is provided)

    # define general and specific integers

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

    # create scorer, if it is empty
    if not scorer:
        seqA_tokens = list(set(seqA))
        seqB_tokens = list(set(seqB))
        for tA in seqA_tokens:
            for tB in seqB_tokens:
                if tA == tB:
                    scorer[tA,tB] = 1
                else:
                    scorer[tA,tB] = -1
    
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



