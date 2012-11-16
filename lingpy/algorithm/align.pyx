"""
Core module for alignment analyses.

"""
__author__ = 'Johann-Mattis List'
__date__ = '2012-16-12'

import random

def _global(
        list listA,
        int lenA,
        list listB,
        int lenB,
        list scorer,
        float scale,
        list almA,
        list almB
        ):
    """
    Internal function for global alignment analyses. 
    """

    cdef int i
    cdef int j
    cdef float gapA
    cdef float gapB
    cdef float match

    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

    traceback[0][0] = 1

    for i in range(1,lenA+1):
        matrix[0][i] = matrix[0][i-1] + scorer[0][i] * scale
        traceback[0][i] = 2
    for i in range(1,lenB+1):
        matrix[i][0] = matrix[i-1][0] + scorer[i][0] * scale
        traceback[i][0] = 3

    for i in range(1,lenB+1):
        for j in range(1,lenA+1):
            if listB[i-1] < 0 and listA[j-1] > 0 and j != lenA:
                gapA = matrix[i-1][j] - 1000000000
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + scorer[i][0] * scale
            else:
                gapA = matrix[i-1][j] + scorer[i][0]
            
            if listA[j-1] < 0 and listB[i-1] > 0 and i != lenB:
                gapB = matrix[i][j-1] - 1000000000
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + scorer[0][j] * scale
            else:
                gapB = matrix[i][j-1] + scorer[0][j]

            match = matrix[i-1][j-1] + scorer[i][j]

            if gapA > match and gapA >= gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    cdef float sim = matrix[lenB][lenA]

    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA[j] = almA[j] + 1
            i -= 1 
        elif traceback[i][j] == 1: 
            i -= 1
            j -= 1
        
        else:
            almB[i] = almB[i] + 1
            j -= 1

    return sim

def _local(
        list listA,
        int lenA,
        list listB,
        int lenB,
        list scorer,
        float scale,
        list almA,
        list almB
        ):
    """
    Internal function for local alignment analyses. 
    """

    cdef int i
    cdef int j
    cdef int k
    cdef float gapA
    cdef float gapB
    cdef float match
    cdef float null
    cdef int imax
    cdef int jmax
    cdef float max_score = 0.0

    cdef list matrix = [[0.0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

    for i in range(1,lenB+1):
        for j in range(1,lenA+1):
            
            null = 0.0

            if listB[i-1] < 0 and listA[j-1] > 0 and j != lenA:
                gapA = matrix[i-1][j] - 1000000000
                null = -1000000000
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + scorer[i][0] * scale
                
            else:
                gapA = matrix[i-1][j] + scorer[i][0]
            
            if listA[j-1] < 0 and listB[i-1] > 0 and i != lenB:
                gapB = matrix[i][j-1] - 1000000000
                null = -1000000000
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + scorer[0][j] * scale
            else:
                gapB = matrix[i][j-1] + scorer[0][j]

            match = matrix[i-1][j-1] + scorer[i][j]

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
                traceback[i][j] = 0

            if matrix[i][j] >= max_score:
                max_score = matrix[i][j]
                imax = i
                jmax = j

    cdef float sim = matrix[imax][jmax]

    i = imax
    j = jmax

    for k in range(j,lenA):
        almA[k] = -1
    for k in range(i,lenB):
        almB[k] = -1

    while traceback[i][j] != 0:
        if traceback[i][j] == 3:
            almA[j] = almA[j] + 1
            i -= 1 
        elif traceback[i][j] == 1: 
            i -= 1
            j -= 1
        
        elif traceback[i][j] == 2:
            almB[i] = almB[i] + 1
            j -= 1
        else:
            break
    for k in range(j):
        almA[k] = -1
    for k in range(i):
        almB[k] = -1

    return sim

def _overlap(
        list listA,
        int lenA,
        list listB,
        int lenB,
        list scorer,
        float scale,
        list almA,
        list almB
        ):
    """
    Internal function for global alignment analyses. 
    """

    cdef int i
    cdef int j
    cdef float gapA
    cdef float gapB
    cdef float match

    cdef list matrix = [[0.0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

    traceback[0][0] = 1

    for i in range(1,lenA+1):
        traceback[0][i] = 2
    for i in range(1,lenB+1):
        traceback[i][0] = 3

    for i in range(1,lenB+1):
        for j in range(1,lenA+1):
            if listB[i-1] < 0 and listA[j-1] > 0 and j != lenA:
                gapA = matrix[i-1][j] - 1000000000
            elif j == lenA:
                gapA = matrix[i-1][j]
            elif traceback[i-1][j] == 3:
                gapA = matrix[i-1][j] + scorer[i][0] * scale
            else:
                gapA = matrix[i-1][j] + scorer[i][0]
            
            if listA[j-1] < 0 and listB[i-1] > 0 and i != lenB:
                gapB = matrix[i][j-1] - 1000000000
            elif i == lenB:
                gapB = matrix[i][j-1]
            elif traceback[i][j-1] == 2:
                gapB = matrix[i][j-1] + scorer[0][j] * scale
            else:
                gapB = matrix[i][j-1] + scorer[0][j]

            match = matrix[i-1][j-1] + scorer[i][j]

            if gapA > match and gapA > gapB:
                matrix[i][j] = gapA
                traceback[i][j] = 3
            elif match >= gapB:
                matrix[i][j] = match
                traceback[i][j] = 1
            else:
                matrix[i][j] = gapB
                traceback[i][j] = 2

    cdef float sim = matrix[lenB][lenA]

    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA[j] = almA[j] + 1
            i -= 1 
        elif traceback[i][j] == 1: 
            i -= 1
            j -= 1
        
        else:
            almB[i] = almB[i] + 1
            j -= 1

    return sim

def _dialign(
        list listA,
        int lenA,
        list listB,
        int lenB,
        list scorer,
        float scale,
        list almA,
        list almB
        ):
    """
    Internal function for global alignment analyses using the DIALIGN algorithm. 
    """

    cdef int i,j,k,l
    cdef int minimum
    cdef float old_score,new_score
    cdef int old_length,new_length

    cdef float scoreA,scoreB,max_score,sim

    cdef list matrix = [[0.0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

    for i in range(1,lenA+1):
        traceback[0][i] = 2
    for i in range(1,lenB+1):
        traceback[i][0] = 3

    for i in range(1,lenB+1):
        for j in range(1,lenA+1):
            
            # determine the minimum sequence
            if i < j:
                minimum = i
            else:
                minimum = j

            old_score = 0.0
            old_length = 1

            for k in range(minimum):
                new_score = matrix[i-k-1][j-k-1]

                for l in range(k,-1,-1):
                    new_score += scorer[i-l][j-l]
                new_length = k+1

                if new_score > old_score:
                    old_score = new_score
                    old_length = new_length

            if listB[i-1] < 0 and listA[j-1] > 0 and j != lenA:
                scoreA = matrix[i-1][j] - 1000000000
            else:
                scoreA = matrix[i-1][j]

            if listA[j-1] < 0 and listB[i-1] > 0 and i != lenB:
                scoreB = matrix[i][j-1] - 1000000000
            else:
                scoreB = matrix[i][j-1]
            
            # determine the maximum score
            if scoreA >= old_score and scoreA > scoreB:
                max_score = scoreA
                traceback[i][j] = 3
            elif old_score > scoreB:
                max_score = old_score
                for k in range(old_length):
                    traceback[i-k][j-k] = 1
            else:
                max_score = scoreB
                traceback[i][j] = 2

            matrix[i][j] = max_score
        
    sim = matrix[lenB][lenA]

    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            almA[j] = almA[j] + 1
            i -= 1 
        elif traceback[i][j] == 1: 
            i -= 1
            j -= 1
        
        else:
            almB[i] = almB[i] + 1
            j -= 1

    return sim

def align_pairwise(
        list seqA,
        list seqB,
        list wghA,
        list wghB,
        list resA,
        list resB,
        str prsA,
        str prsB,
        dict score_dict,
        float scale,
        float sonority_factor,
        str mode
        ):
    """
    Basic function for alignment analyses. 
    """
    cdef int lA = len(seqA)
    cdef int lB = len(seqB)

    cdef list outA = seqA[:]
    cdef list outB = seqB[:]

    cdef int k
    cdef int l

    if mode == "global":
        aligner = _global
    elif mode == "local":
        aligner = _local
    elif mode == "overlap":
        aligner = _overlap
    elif mode == "dialign":
        aligner = _dialign
    
    cdef list scorer = [[0 for i in range(lA+1)] for j in range(lB+1)]

    for k in range(1,lB+1):
        scorer[k][0] = wghB[k-1]
    for k in range(1,lA+1):
        scorer[0][k] = wghA[k-1]

    for k in range(1,lB+1):
        for l in range(1,lA+1):
            score = score_dict[seqA[l-1],seqB[k-1]]

            if prsA[l-1] == prsB[k-1]: # and prsA[l-1].lower() == prsB[k-1].lower():
                score = score * (1.0 + sonority_factor)
            scorer[k][l] = score

    cdef list listA = [k for k in resA]
    cdef list listB = [k for k in resB]

    cdef list almA = [0 for k in range(lA+1)]
    cdef list almB = [0 for k in range(lB+1)]

    cdef float sim = aligner(
            listA,
            lA,
            listB,
            lB,
            scorer,
            scale,
            almA,
            almB
            )

    for k in range(lA,-1,-1):
        if almA[k] > 0:
            for l in range(almA[k]):
                outA.insert(k,"-")
        elif almA[k] < 0:
            outA[k] = "*"
    for k in range(lB,-1,-1):
        if almB[k] > 0:
            for l in range(almB[k]):
                outB.insert(k,"-")
        elif almB[k] < 0:
            outB[k] = "*"

    return outA,outB,sim

def align_sequences_pairwise(
        list seqs,
        list weights,
        list restrictions,
        list prosodics,
        dict score_dict,
        float scale,
        float sonority_factor,
        str mode
        ):
    """
    Function takes a list of sequences as input and returns all possible
    pairwise alignments between all sequences.
    """

    cdef int lS = len(seqs)

    cdef int i,j,k,l

    cdef float score
    cdef float sim

    cdef list alignments = []

    # more and more cdefs...
    cdef list seqA,seqB,wghA,wghB,resA,resB,outA,outB,almA,almB
    #cdef str prsA,prsB
    cdef int lA,lB
    cdef list scorer,listA,listB

    if mode == "global":
        aligner = _global
    elif mode == "local":
        aligner = _local
    elif mode == "overlap":
        aligner = _overlap
    elif mode == "dialign":
        aligner = _dialign

    for i in range(lS):
        for j in range(lS):
            if i < j:
                seqA = seqs[i]
                seqB = seqs[j]
                wghA = weights[i]
                wghB = weights[j]
                resA = restrictions[i]
                resB = restrictions[j]
                prsA = prosodics[i]
                prsB = prosodics[j]

                lA = len(seqA)
                lB = len(seqB)

                outA = seqA[:]
                outB = seqB[:]

                scorer = [[0.0 for k in range(lA+1)] for l in range(lB+1)]
                
                for k in range(1,lB+1):
                    scorer[k][0] = wghB[k-1]
                for k in range(1,lA+1):
                    scorer[0][k] = wghA[k-1]

                for k in range(1,lB+1):
                    for l in range(1,lA+1):
                        score = score_dict[seqA[l-1],seqB[k-1]]

                        if prsA[l-1] == prsB[k-1]: # and prsA[l-1].lower() == prsB[k-1].lower():
                            score = score * (1.0 + sonority_factor)
                        scorer[k][l] = score

                listA = [k for k in resA]
                listB = [k for k in resB]
                almA = [0 for k in range(lA+1)]
                almB = [0 for k in range(lB+1)]

                sim = aligner(
                        listA,
                        lA,
                        listB,
                        lB,
                        scorer,
                        scale,
                        almA,
                        almB
                        )

                for k in range(lA,-1,-1):
                    if almA[k] > 0:
                        for l in range(almA[k]):
                            outA.insert(k,"-")
                    elif almA[k] < 0:
                        outA[k] = "*"

                for k in range(lB,-1,-1):
                    if almB[k] > 0:
                        for l in range(almB[k]):
                            outB.insert(k,"-")
                    elif almB[k] < 0:
                        outB[k] = "*"

                alignments.append((outA,outB,sim))

    return alignments

def align_sequence_pairs(
        list seqs,
        list weights,
        list restrictions,
        list prosodics,
        dict score_dict,
        float scale,
        float sonority_factor,
        str mode
        ):
    """
    Function takes a list of sequence pairs as input and returns the aligned
    sequence pairs.
    """

    cdef int lS = len(seqs)

    cdef int i,j,k,l

    cdef float score
    cdef float sim

    cdef list alignments = []

    # more and more cdefs...
    cdef list seqA,seqB,wghA,wghB,resA,resB,outA,outB,almA,almB
    # cdef str prsA,prsB
    cdef int lA,lB
    cdef list scorer,listA,listB

    if mode == "global":
        aligner = _global
    elif mode == "local":
        aligner = _local
    elif mode == "overlap":
        aligner = _overlap
    elif mode == "dialign":
        aligner = _dialign

    for i in range(lS):
        seqA = seqs[i][0]
        seqB = seqs[i][1]
        wghA = weights[i][0]
        wghB = weights[i][1]
        resA = restrictions[i][0]
        resB = restrictions[i][1]
        prsA = prosodics[i][0]
        prsB = prosodics[i][1]

        lA = len(seqA)
        lB = len(seqB)

        outA = seqA[:]
        outB = seqB[:]

        scorer = [[0.0 for k in range(lA+1)] for l in range(lB+1)]

        for k in range(1,lB+1):
            scorer[k][0] = wghB[k-1]
        for k in range(1,lA+1):
            scorer[0][k] = wghA[k-1]

        for k in range(1,lB+1):
            for l in range(1,lA+1):
                score = score_dict[seqA[l-1],seqB[k-1]]

                if prsA[l-1] == prsB[k-1]: # and prsA[l-1].lower() == prsB[k-1].lower():
                    score = score * (1.0 + sonority_factor)
                scorer[k][l] = score

        listA = [k for k in resA]
        listB = [k for k in resB]
        almA = [0 for k in range(lA+1)]
        almB = [0 for k in range(lB+1)]

        sim = aligner(
                listA,
                lA,
                listB,
                lB,
                scorer,
                scale,
                almA,
                almB
                )
        for k in range(lA,-1,-1):
            if almA[k] > 0:
                for l in range(almA[k]):
                    outA.insert(k,"-")
            elif almA[k] < 0:
                outA[k] = "*"

        for k in range(lB,-1,-1):
            if almB[k] > 0:
                for l in range(almB[k]):
                    outB.insert(k,"-")
            elif almB[k] < 0:
                outB[k] = "*"

        alignments.append((outA,outB,sim))

    return alignments

def random_align_sequence_pairs(
        list seqs,
        list weights,
        list restrictions,
        list prosodics,
        dict score_dict,
        float scale,
        float sonority_factor,
        str mode,
        int runs
        ):
    """
    Function takes a list of sequences pairs as input and returns a dictionary
    of correspondence frequencies.
    """

    cdef int lS = len(seqs)

    cdef int i,j,k,l,n,run

    cdef float score
    cdef float sim

    # more and more cdefs...
    cdef list seqA,seqB,wghA,wghB,resA,resB,outA,outB,almA,almB
    #cdef str prsA,prsB
    cdef int lA,lB
    cdef list scorer,listA,listB
    cdef dict corrs = {}

    if mode == "global":
        aligner = _global
    elif mode == "local":
        aligner = _local
    elif mode == "overlap":
        aligner = _overlap
    elif mode == "dialign":
        aligner = _dialign

    cdef list randoms = [i for i in range(lS)]
    cdef dict alm_pairs = {}

    for i in range(lS):
        for run in range(runs):
            random.shuffle(randoms)
            
            j = randoms[i]

            if (i,j) not in alm_pairs:
                seqA = seqs[i][0]
                seqB = seqs[j][1]
                wghA = weights[i][0]
                wghB = weights[j][1]
                resA = restrictions[i][0]
                resB = restrictions[j][1]
                prsA = prosodics[i][0]
                prsB = prosodics[j][1]

                lA = len(seqA)
                lB = len(seqB)

                outA = seqA[:]
                outB = seqB[:]

                scorer = [[0.0 for k in range(lA+1)] for l in range(lB+1)]

                for k in range(1,lB+1):
                    scorer[k][0] = wghB[k-1]
                for k in range(1,lA+1):
                    scorer[0][k] = wghA[k-1]

                for k in range(1,lB+1):
                    for l in range(1,lA+1):
                        score = score_dict[seqA[l-1],seqB[k-1]]

                        if prsA[l-1] == prsB[k-1]: # and prsA[l-1].lower() == prsB[k-1].lower():
                            score = score * (1.0 + sonority_factor)
                        scorer[k][l] = score

                listA = [k for k in resA]
                listB = [k for k in resB]
                almA = [0 for k in range(lA+1)]
                almB = [0 for k in range(lB+1)]

                sim = aligner(
                        listA,
                        lA,
                        listB,
                        lB,
                        scorer,
                        scale,
                        almA,
                        almB
                        )
                for k in range(lA,-1,-1):
                    if almA[k] > 0:
                        for l in range(almA[k]):
                            outA.insert(k,"-")
                    elif almA[k] < 0:
                        outA[k] = "*"

                for k in range(lB,-1,-1):
                    if almB[k] > 0:
                        for l in range(almB[k]):
                            outB.insert(k,"-")
                    elif almB[k] < 0:
                        outB[k] = "*"

                alm_pairs[i,j] = (outA,outB)
            else:
                outA = alm_pairs[i,j][0]
                outB = alm_pairs[i,j][1]
            
            if mode == 'local':
                outA = [m for m in outA if m != '*']
                outB = [m for m in outB if m != '*']

            for n in range(len(outA)):
                if (outA[n],outB[n]) in corrs:
                    corrs[outA[n],outB[n]] += 1.0 / runs
                else:
                    corrs[outA[n],outB[n]] = 1.0 / runs

    return corrs

def edit_dist(
        list seqA,
        list seqB
        ):
    """
    Return the normalized edit-distance between two strings.
    """
    
    cdef int lenA = len(seqA)
    cdef int lenB = len(seqB)
    cdef int gapA,gapB,match
    cdef int i,j
    
    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    
    for i in range(1,lenA+1):
        matrix[0][i] = i
    for i in range(1,lenB+1):
        matrix[i][0] = i

    for i in range(1,lenB+1):
        for j in range(1,lenA+1):
            
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

    cdef float sim = matrix[lenB][lenA]
    cdef float dist = sim / float(max([lenA,lenB]))

    return dist

def nw_align(
        list seqA,
        list seqB,
        dict scorer,
        int gap = -1
        ):
    """
    Align two sequences using the Needleman-Wunsch algorithm.
    """
    
    # get the lengths of the strings
    cdef int lenA = len(seqA)
    cdef int lenB = len(seqB)

    # define lists for tokens (in case no scoring function is provided)
    cdef list seqA_tokens,seqB_tokens
    cdef str tA,tB

    # define general and specific integers
    cdef int i,j
    cdef int sim # stores the similarity score

    # define values for the main loop
    cdef int gapA,gapB,match,penalty # for the loop
 
    # define values for the traceback
    cdef list almA = seqA[:]
    cdef list almB = seqB[:] 
    cdef str gap_char = '-' # the gap character

    # create matrix and traceback
    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

    # initialize matrix and traceback
    for i in range(1,lenA+1):
        matrix[0][i] = matrix[0][i-1] + gap
        traceback[0][i] = 2
    for i in range(1,lenB+1):
        matrix[i][0] = matrix[i-1][0] + gap
        traceback[i][0] = 3

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
            almA.insert(j,gap_char)
            i -= 1
        elif traceback[i][j] == 1:
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:
            almB.insert(i,gap_char)
            j -= 1
        else:
            break

    # return the alignment as a tuple of prefix, alignment, and suffix
    return (almA,almB,sim)


def sw_align(
        list seqA,
        list seqB,
        dict scorer,
        int gap = -1
        ):
    """
    Align two sequences using the Smith-Waterman algorithm.
    """
    
    # get the lengths of the strings
    cdef int lenA = len(seqA)
    cdef int lenB = len(seqB)

    # define lists for tokens (in case no scoring function is provided)
    cdef list seqA_tokens,seqB_tokens
    cdef str tA,tB

    # define general and specific integers
    cdef int i,j
    cdef int sim # stores the similarity score

    # define values for the main loop
    cdef int gapA,gapB,match,penalty # for the loop
    cdef int null = 0 # constant during the loop
    cdef int imax = 1 # for the loop
    cdef int jmax = 1 # for the loop
    cdef int max_score = 0 # for the loo

    # define values for the traceback
    cdef int igap = 0
    cdef int jgap = 0 
    cdef list almA = seqA[:]
    cdef list almB = seqB[:] 
    cdef str gap_char = '-' # the gap character

    # create matrix and traceback
    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

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

    # return the alignment as a tuple of prefix, alignment, and suffix
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
        list seqA,
        list seqB,
        dict scorer,
        int gap = -1
        ):
    """
    Align two sequences using the Waterman-Eggert algorithm.
    """
    
    # get the lengths of the strings
    cdef int lenA = len(seqA)
    cdef int lenB = len(seqB)

    # define lists for tokens (in case no scoring function is provided)
    cdef list seqA_tokens,seqB_tokens
    cdef str tA,tB

    # define general and specific integers
    cdef int i,j
    cdef int sim # stores the similarity score

    # define values for the main loop
    cdef int gapA,gapB,match,penalty # for the loop
    cdef int null = 0 # constant during the loop
    cdef int imax,jmax # for the loop
    cdef int imin,jmin
    cdef int max_score # for the loo

    # define values for the traceback
    cdef int igap = 0
    cdef int jgap = 0 
    cdef list almA,almB 
    cdef str gap_char = '-' # the gap character

    # create a tracer for positions in the matrix
    cdef list tracer = [0 for i in range(lenA+1)]
    cdef int idx

    # create matrix and traceback
    cdef list matrix = [[0 for i in range(lenA+1)] for j in range(lenB+1)]
    cdef list traceback = [[0 for i in range(lenA+1)] for j in range(lenB+1)]

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

    
    # make list of alignments
    cdef list out = []

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
        almA = seqA[:]
        almB = seqB[:]

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

    # return the alignment as a tuple of prefix, alignment, and suffix
    return out 

