# we start with basic alignment functions
def globalign(
        list seqA,
        list seqB,
        int M, # length of seqA
        int N, # length of seqB
        int gop,
        float scale,
        object scorer
        ):
    """
    Carry out global alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as lists.
    M, N : int
        The length of the two sequences.
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.

    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity score.
    
    Notes
    -----
    This algorithm carries out classical Needleman-Wunsch alignment
    (:evobib:`Needleman1970`).
    
    See also
    --------
    ~lingpy.algorithm.cython.talign.semi_globalign
    ~lingpy.algorithm.cython.talign.localign
    ~lingpy.algorithm.cython.talign.dialign

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
        matrix[0][i] = matrix[0][i-1] + gop * scale
        traceback[0][i] = 2
    for i in range(1,N+1):
        matrix[i][0] = matrix[i-1][0] + gop * scale
        traceback[i][0] = 3

    # start the loop
    for i in range(1,N+1):
        for j in range(1,M+1):

            # calculate costs for gapA
            if traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + gop * scale
            else:
                gapA = matrix[i-1][j] + gop

            # calculate costs for gapB
            if traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gop * scale
            else:
                gapB = matrix[i][j-1] + gop

            # get the score
            match = matrix[i-1][j-1] + scorer[seqA[j-1],seqB[i-1]]

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
        int M, # length of seqA
        int N, # length of seqB
        int gop,
        float scale,
        object scorer
        ):
    """
    Carry out semi-global alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as lists.
    M, N : int
        The length of the two sequences.
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.
    
    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity score.

    Notes
    -----
    This algorithm carries out semi-global alignment 
    (:evobib:`Durbin2002`).
    
    See also
    --------
    ~lingpy.algorithm.cython.talign.globalign
    ~lingpy.algorithm.cython.talign.localign
    ~lingpy.algorithm.cython.talign.dialign

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
                gapA = matrix[i-1][j] + gop * scale
            else:
                gapA = matrix[i-1][j] + gop

            # calculate costs for gapB
            if i == N:
                gapB = matrix[i][j-1]
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gop * scale
            else:
                gapB = matrix[i][j-1] + gop

            # calculate costs for match

            # get the score
            match = matrix[i-1][j-1] + scorer[seqA[j-1],seqB[i-1]]

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
        int M, # length of seqA
        int N, # length of seqB
        int gop,
        float scale,
        object scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    
    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as lists.
    M, N : int
        The length of the two sequences.
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.
    
    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity score.

    Notes
    -----
    This algorithm carries out local alignment 
    (:evobib:`Smith1981`).
    
    See also
    --------
    ~lingpy.algorithm.cython.talign.globalign
    ~lingpy.algorithm.cython.talign.semi_globalign
    ~lingpy.algorithm.cython.talign.dialign
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
                gapA = matrix[i-1][j] + gop * scale
            else:
                gapA = matrix[i-1][j] + gop

            # calculate costs for gapB
            if traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + gop * scale
            else:
                gapB = matrix[i][j-1] + gop

            # calculate costs for match

            # get the score
            match = matrix[i-1][j-1] + scorer[seqA[j-1],seqB[i-1]]

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
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        object scorer
        ):
    """
    Carry out dialign alignment of two sequences.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as lists.
    M, N : int
        The length of the two sequences.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.
    
    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity score.

    Notes
    -----
    This algorithm carries out dialign alignment 
    (:evobib:`Morgenstern1996`).
    
    See also
    --------
    ~lingpy.algorithm.cython.talign.globalign
    ~lingpy.algorithm.cython.talign.semi_globalign
    ~lingpy.algorithm.cython.talign.localign
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
                    match += scorer[seqA[j-l-1],seqB[i-l-1]]

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
        int gop,
        float scale,
        object scorer,
        str mode,
        int distance = 0
        ):
    """
    Align a pair of sequences.

    Parameters
    ----------
    seqA, seqB : list
        The sequences to be aligned, passed as lists.
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.
    mode : { "global", "local", "overlap", "dialign" }
        Select the mode for the alignment analysis ("overlap" refers to
        semi-global alignments).
    distance : int (default=0)
        Select whether you want distances or similarities to be returned (0
        indicates similarities, 1 indicates distances, 2 indicates both).

    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity score.

    Notes
    -----
    This is a utility function that allows calls any of the four classical
    alignment functions (:py:class:`lingpy.algorithm.cython.talign.globalign`
    :py:class:`lingpy.algorithm.cython.talign.semi_globalign`,
    :py:class:`lingpy.algorithm.cython.talign.lotalign`,
    :py:class:`lingpy.algorithm.cython.talign.dialign`,) and their secondary counterparts.

    See also
    --------
    ~lingpy.algorithm.cython.talign.align_pairwise
    ~lingpy.algorithm.cython.talign.align_pairs

    Returns
    -------
    alignment : tuple
        The aligned sequences and the similarity or distance scores, or both.

    """
    # define basic types
    cdef int i
    cdef list almA,almB
    cdef float sim,dist,simA,simB

    # get length of seqA,seqB
    cdef int M = len(seqA)
    cdef int N = len(seqB)

    # determine the mode
    if mode == "global":
        
        # carry out the alignment
        almA,almB,sim = globalign(
                seqA,
                seqB,
                M,
                N,
                gop,
                scale,
                scorer
                )

    elif mode == "local":
        
        # carry out the alignment
        almA,almB,sim = localign(
                seqA,
                seqB,
                M,
                N,
                gop,
                scale,
                scorer
                )

    elif mode == "overlap":
        
        # carry out the alignment
        almA,almB,sim = semi_globalign(
                seqA,
                seqB,
                M,
                N,
                gop,
                scale,
                scorer
                )

    elif mode == "dialign":
        almA,almB,sim = dialign(
                seqA,
                seqB,
                M,
                N,
                scale,
                scorer
                )

    # calculate distance, if this is needed
    if distance > 0:
        simA = sum([scorer[seqA[i],seqA[i]] for i in range(M)])
        simB = sum([scorer[seqB[i],seqB[i]] for i in range(N)])

        dist = 1 - ( ( 2 * sim ) / ( simA + simB ) )
        if distance == 1:
            return almA,almB,dist
        else:
            return almA,almB,sim,dist
    else:
        return almA,almB,sim
    
def align_pairwise(
        list seqs,
        int gop,
        float scale,
        object scorer,
        str mode
        ):
    """
    Align all sequences pairwise.

    Parameters
    ----------
    seqs : list
        The sequences to be aligned, passed as lists.
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.
    mode : { "global", "local", "overlap", "dialign" }
        Select the mode for the alignment analysis ("overlap" refers to
        semi-global alignments).

    Returns
    -------
    alignments : list
        A list of tuples, containing the aligned sequences, the similarity
        and the distance scores.

    Notes
    -----
    This function aligns all possible pairs between the sequences you pass to
    it. It is important for multiple alignment, where it can be used to
    construct the guide tree.

    See also
    --------
    ~lingpy.algorithm.cython.talign.align_pair
    ~lingpy.algorithm.cython.talign.align_pairs
    """
    # define basic stuff
    cdef list alignments = []
    cdef int lS = len(seqs)
    
    cdef int i,j,k,lenA,lenB
    cdef list almA,almB,seqA,seqB
    cdef float sim,simA,simB,dist

    # get self-scores
    cdef list sims = [0.0 for i in range(lS)]
    cdef list lens = [0 for i in range(lS)]

    for i in range(lS):
        seqA = seqs[i]
        k = len(seqA)
        sim = sum([scorer[seqA[j],seqA[j]] for j in range(k)])
        lens[i] = k
        sims[i] = sim
    
    if mode == "global":
        # start loop
        for i in range(lS):
            for j in range(lS):
                if i < j:
                    seqA,seqB = seqs[i],seqs[j]
                    simA,simB = sims[i],sims[j]
                    lenA,lenB = lens[i],lens[j]
                    almA,almB,sim = globalign(
                            seqA,
                            seqB,
                            lenA,
                            lenB,
                            gop,
                            scale,
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
                    simA,simB = sims[i],sims[j]
                    lenA,lenB = lens[i],lens[j]

                    # check for secondary structures
                    almA,almB,sim = localign(
                            seqA,
                            seqB,
                            lenA,
                            lenB,
                            gop,
                            scale,
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
                    simA,simB = sims[i],sims[j]
                    lenA,lenB = lens[i],lens[j]

                    almA,almB,sim = semi_globalign(
                            seqA,
                            seqB,
                            lenA,
                            lenB,
                            gop,
                            scale,
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
                    simA,simB = sims[i],sims[j]
                    lenA,lenB = lens[i],lens[j]

                    almA,almB,sim = dialign(
                           seqA,
                           seqB,
                           lenA,
                           lenB,
                           scale,
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

    return alignments
    
def align_pairs(
        list seqs,
        int gop,
        float scale,
        object scorer,
        str mode,
        int distance = 0
        ):
    """
    Align multiple sequence pairs.

    Parameters
    ----------
    seqs : list
        The sequences to be aligned, passed as lists.
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale.
    scorer : { dict, ~lingpy.algorithm.cython.misc.ScoreDict }
        The scoring dictionary containing scores for all possible segment
        combinations in the two sequences.
    mode : { "global", "local", "overlap", "dialign" }
        Select the mode for the alignment analysis ("overlap" refers to
        semi-global alignments).
    distance : int (default=0)
        Indicate whether distances or similarities should be returned.

    Returns
    -------
    alignments : list
        A list of tuples, containing the aligned sequences, and the similarity
        or the distance scores.

    Notes
    -----
    This function aligns all pairs which are passed to
    it. 

    See also
    --------
    ~lingpy.algorithm.cython.talign.align_pair
    ~lingpy.algorithm.cython.talign.align_pairwise

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

        if mode == "global":
            almA,almB,sim = globalign(
                   seqA,
                   seqB,
                   M,
                   N,
                   gop,
                   scale,
                   scorer
                   )
        elif mode == "local":
            almA,almB,sim = localign(
                   seqA,
                   seqB,
                   M,
                   N,
                   gop,
                   scale,
                   scorer
                   )

        elif mode == "overlap":
            almA,almB,sim = semi_globalign(
                   seqA,
                   seqB,
                   M,
                   N,
                   gop,
                   scale,
                   scorer
                   )

        elif mode == "dialign":
            almA,almB,sim = dialign(
                   seqA,
                   seqB,
                   M,
                   N,
                   scale,
                   scorer
                   )

        # calculate distances if option is chose
        if distance > 0:
            simA = sum([scorer[seqA[i],seqA[i]] for i in range(M)])
            simB = sum([scorer[seqB[i],seqB[i]] for i in range(N)])

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
        int gop,
        float scale,
        object scorer,
        str mode,
        float gap_weight
        ):
    """
    Align two profiles using the basic modes.

    Parameters
    ----------
    profileA, profileB : list
        Two-dimensional list for each of the profiles. 
    gop : int
        The gap opening penalty.
    scale : float
        The gap extension scale by which consecutive gaps are reduced. LingPy
        uses a scale rather than a constant gap extension penalty. 
    scorer : { dict, :py:class:`lingpy.algorithm.cython.misc.ScoreDict` }
        The scoring function which needs to provide scores for all
        segments in the two profiles.
    mode : { "global", "overlap", "dialign" }
        Select one of the four basic modes for alignment analyses.
    gap_weight : float
        This handles the weight that is given to gaps in a column. If you set
        it to 0, for example, this means that all gaps will be ignored when
        determining the score for two columns in the profile.

    Notes
    -----
    This function computes alignments of two profiles of multiple sequences
    (see :evobib:`Durbin2002` for details on profiles)
    and is important for multiple alignment analyses.

    Returns
    -------
    alignment : tuple
        The aligned profiles, and the overall similarity of the profiles.

    See also
    --------
    ~lingpy.algorithm.cython.talign.score_profile
    ~lingpy.algorithm.cython.talign.swap_score_profile
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
    
    if mode == "global":
        almA,almB,sim = globalign(
                listA,
                listB,
                M,
                N,
                gop,
                scale,
                tmp_scorer
                )
    elif mode == "overlap":
        almA,almB,sim = semi_globalign(
                listA,
                listB,
                M,
                N,
                gop,
                scale,
                tmp_scorer
                )
    elif mode == "dialign":
        almA,almB,sim = dialign(
                listA,
                listB,
                M,
                N,
                scale,
                tmp_scorer
                )
     
    return almA,almB,sim

def score_profile(
        list colA,
        list colB,
        object scorer,
        int gop,
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
    ~lingpy.algorithm.cython.talign.align_profile
    ~lingpy.algorithm.cython.talign.swap_score_profile
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
            elif charA == 'X' and charB == 'X':
                counter += gap_weight
            else:
                score += float(gop)
                counter += 1.0
    return score / counter

def swap_score_profile(
        list colA,
        list colB,
        object scorer,
        float gap_weight = 0.0,
        int swap_penalty = -1
        ):
    """
    Basic function for the scoring of profiles in swapped sequences.

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
    ~lingpy.algorithm.cython.talign.align_profile
    ~lingpy.algorithm.cython.talign.score_profile

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


