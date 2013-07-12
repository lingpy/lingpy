# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-10 18:33
# modified : 2013-03-10 18:35
"""
Basic module for sound-class based alignment analyses.
"""

__author__="Johann-Mattis List"
__date__="2013-03-10"

# we start with basic alignment functions
def globalign(
        seqA,
        seqB,
        M, # length of seqA
        N, # length of seqB
        gop,
        scale,
        scorer
        ):
    """
    Carry out global alignment of two sequences.
    """

    # declare integers
#     cdef int i,j

    # declare floats
#     cdef gapA,gapB,match,sim

    # declare lists
    almA = []
    almB = []

    # create matrix and traceback
    matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]

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
        seqA,
        seqB,
        M, # length of seqA
        N, # length of seqB
        gop,
        scale,
        scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
#     cdef int i,j

    # declare floats
#     cdef gapA,gapB,match,sim

    # declare lists
    almA = []
    almB = []

    # create matrix and traceback
    matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]

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
        seqA,
        seqB,
        M, # length of seqA
        N, # length of seqB
        gop,
        scale,
        scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
#     cdef int i,j,k,l

    # declare floats
#     cdef gapA,gapB,match,sim

    # declare char-character
#     cdef str x

    # declare lists
    almA = []
    almB = []

    # create matrix and traceback
    matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]

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
        seqA,
        seqB,
        M, # length of seqA
        N, # length of seqB
        scale,
        scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
#     cdef int i,j,k,l,o,p

    # declare floats
#     cdef gapA,gapB,match,sim,tmp_match

    # declare lists
    almA = []
    almB = []

    # create matrix and traceback
    matrix = [[0.0 for i in range(M+1)] for j in range(N+1)]
    traceback = [[0 for i in range(M+1)] for j in range(N+1)]

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
        seqA,
        seqB,
        gop,
        scale,
        scorer,
        mode,
        distance = 0
        ):
    """
    Align a pair of sequences.
    """
    # define basic types
#     cdef int i
#     cdef list almA,almB
#     cdef float sim,dist,simA,simB

    # get length of seqA,seqB
    M = len(seqA)
    N = len(seqB)

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
        seqs,
        gop,
        scale,
        scorer,
        mode
        ):
    """
    Align a list of sequences pairwise.
    """
    # define basic stuff
    alignments = []
    lS = len(seqs)
    
#     cdef int i,j,k,lenA,lenB
#     cdef list almA,almB,seqA,seqB
#     cdef float sim,simA,simB,dist

    # get self-scores
    sims = [0.0 for i in range(lS)]
    lens = [0 for i in range(lS)]

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
        seqs,
        gop,
        scale,
        scorer,
        mode,
        distance = 0
        ):
    """
    Align multiple sequence pairs.
    """
    # basic defs
#     cdef int i,j,M,N,lP
#     cdef list seqA,seqB,almA,almB
#     cdef float sim
    alignments = []

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
        profileA,
        profileB,
        gop,
        scale,
        scorer,
        mode,
        gap_weight
        ):
    """
    Align two profiles using the basic modes.
    """

    # basic defs
#     cdef int i,j,k,l,M,N,O,P
#     cdef float sim,count
#     cdef str charA,charB
#     cdef list listA,listB,almA,almB
    
    M = len(profileA)
    N = len(profileB)
    O = len(profileA[0])
    P = len(profileB[0])

    tmp_scorer = {}

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
                gop,
                scale,
                tmp_scorer
                )
     
    return almA,almB,sim

def score_profile(
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
#     cdef int i,j
#     cdef str charA,charB

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

def swap_score_profile(
        colA,
        colB,
        scorer,
        gap_weight = 0.0,
        swap_penalty = -1
        ):
    """
    Basic function for the scoring of profiles.
    """
    # basic definitions
#     cdef int i,j
#     cdef str charA,charB

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


