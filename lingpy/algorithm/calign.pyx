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
        dict scorer
        ):
    """
    Carry out global alignment of two sequences.
    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef gapA,gapB,match,sim

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
        traceback[i][0] = 2

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
        dict scorer,
        str r # restricted_chars
        ):
    """
    Carry out global alignment of two sequences.
    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef gapA,gapB,match,sim

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
        traceback[i][0] = 2

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
        dict scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef gapA,gapB,match,sim

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
        traceback[i][0] = 2

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

def secondary_dialign(
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
        dict scorer,
        str r # restricted_chars
        ):
    """
    Carry out global alignment of two sequences.
    """

    # declare integers
    cdef int i,j

    # declare floats
    cdef gapA,gapB,match,sim

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
        traceback[i][0] = 2

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
        dict scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
    cdef int i,j,k,l

    # declare floats
    cdef gapA,gapB,match,sim

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
    return almA,almB,sim

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
        dict scorer,
        str r # restricted_chars
        ):
    """
    Carry out global alignment of two sequences.
    """

    # declare integers
    cdef int i,j,k,l

    # declare floats
    cdef gapA,gapB,match,sim

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
                gapA = matrix[i][j-1] - 1000000
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
    return almA,almB,sim

def dialign(
        list seqA,
        list seqB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        dict scorer
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
    cdef int i,j,k,l,o,p

    # declare floats
    cdef gapA,gapB,match,sim,tmp_match

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
        traceback[i][0] = 2

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
                    tmp_match = scorer[seqA[j-l],seqB[i-l]]

                    # check for common prostrings
                    if proA[j-l] == proB[i-l]:
                        tmp_match += tmp_match * factor
                    elif abs(ord(proA[j-l]) - ord(proB[i-l])) <= 2:
                        tmp_match += tmp_match * factor / 2
                    
                    # get match
                    match = tmp_match

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

def secondary_semi_globalign(
        list seqA,
        list seqB,
        str proA,
        str proB,
        int M, # length of seqA
        int N, # length of seqB
        float scale,
        float factor,
        dict scorer,
        str r # restricted chars
        ):
    """
    Carry out semi-global alignment of two sequences.
    """

    # declare integers
    cdef int i,j,k,l,o,p

    # declare floats
    cdef gapA,gapB,match,sim,tmp_match

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
        traceback[i][0] = 2

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
                    tmp_match = scorer[seqA[j-l],seqB[i-l]]

                    # check for common prostrings
                    if proA[j-l] == proB[i-l]:
                        tmp_match += tmp_match * factor
                    elif proA[j-l] in r and proB[i-l] not in r:
                        tmp_match += -1000000
                    elif proA[j-l] not in r and proB[i-l] in r:
                        tmp_match += -1000000
                    elif abs(ord(proA[j-l]) - ord(proB[i-l])) <= 2:
                        tmp_match += tmp_match * factor / 2
                    
                    # get match
                    match = tmp_match

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
        dict scorer,
        str mode,
        str restricted_chars,
        int distance = 0
        ):
    """
    Align a pair of sequences.
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
        str pros,
        int gop,
        float scale,
        float factor,
        dict scorer,
        str restricted_chars,
        str mode
        ):
    """
    Align a list of sequences pairwise.
    """
    # define basic stuff
    cdef list alignments = []
    cdef int lS = len(seqs)
    
    cdef int i,j,k,lenA,lenB
    cdef list almA,almB,seqA,seqB,gopA,gopB
    cdef float sim,simA,simB,dist
    cdef str proA,proB

    # get self-scores
    cdef list sims = []
    cdef list lens
    for i in range(lS):
        seqA = seqs[i]
        k = len(seqA)
        sim = sum([(1 + factor) * scorer[seqA[j],seqA[j]] for j in range(k)])
        lens += [k]
        sims += [sim]
        gops[i] = [gop * gops[j] for j in range(k)]
    
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

                    # check for secondary structures
                    if not set(restricted_chars).intersection(set(proA+proB)):

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

                    else:
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

                    # check for secondary structures
                    if not set(restricted_chars).intersection(set(proA+proB)):
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

                    else:
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

                    # check for secondary structures
                    if not set(restricted_chars).intersection(set(proA+proB)):
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

                    else:
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

                    # check for secondary structures
                    if not set(restricted_chars).intersection(set(proA+proB)):
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

                    else:
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
    



