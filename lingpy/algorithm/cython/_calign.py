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
        gopA,
        gopB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
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
        seqA,
        seqB,
        gopA,
        gopB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
        scorer,
        r # restricted_chars
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
        seqA,
        seqB,
        gopA,
        gopB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
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
        seqA,
        seqB,
        gopA,
        gopB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
        scorer,
        r # restricted_chars
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
        seqA,
        seqB,
        gopA,
        gopB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
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
        seqA,
        seqB,
        gopA,
        gopB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
        scorer,
        r # restricted_chars
        ):
    """
    Carry out global alignment of two sequences.
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
        seqA,
        seqB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
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
        seqA,
        seqB,
        proA,
        proB,
        M, # length of seqA
        N, # length of seqB
        scale,
        factor,
        scorer,
        r # restricted chars
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
        mode,
        restricted_chars,
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
        seqs,
        gops,
        pros,
        gop,
        scale,
        factor,
        scorer,
        restricted_chars,
        mode
        ):
    """
    Align a list of sequences pairwise.
    """
    # define basic stuff
    alignments = []
    lS = len(seqs)
    
#     cdef int i,j,k,lenA,lenB
#     cdef list almA,almB,seqA,seqB,gopA,gopB
#     cdef float sim,simA,simB,dist
#     cdef str proA,proB

    # get self-scores
    sims = [0.0 for i in range(lS)]
    lens = [0 for i in range(lS)]

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
        seqs,
        gops,
        pros,
        gop,
        scale,
        factor,
        scorer,
        mode,
        restricted_chars,
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
        profileA,
        profileB,
        gopA,
        gopB,
        proA,
        proB,
        gop,
        scale,
        factor,
        scorer,
        restricted_chars,
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
        colA,
        colB,
        scorer,
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
            else:
                counter += gap_weight

    return score / counter

def swap_score_profile(
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

def corrdist(
        threshold,
        seqs,
        gops,
        pros,
        gop,
        scale,
        factor,
        scorer,
        mode,
        restricted_chars
        ):
    """
    Create a correspondence distribution for a given language pair.
    """

    # basic defs
#     cdef int i,j,M,N,lP,l
#     cdef list seqA,seqB,almA,almB
#     cdef float sim
    corrs = {}

    # return number of sequences considered for initial distribution
    included = 0

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

    return corrs,included


