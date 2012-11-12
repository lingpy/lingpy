"""
Core module for alignment analyses.

"""
__author__ = 'Johann-Mattis List'
__date__ = '2012-11-12'

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

