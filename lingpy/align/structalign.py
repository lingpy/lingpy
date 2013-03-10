# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-10 22:28
# modified : 2013-03-10 22:28
"""
Structural Alignment Algorithm
"""

__author__="Johann-Mattis List"
__date__="2013-03-10"

from sys import argv

def structalign(seqA,seqB):
    """
    Carry out a structural alignment analysis using Dijkstra's algorithm.
    """
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

        if not restA and not restB:
            out += [(''.join([a[0] for a in alm]),''.join([a[1] for a in alm]))]

        # start adding match
        if restA and restB:
            residues = (restA[0],restB[0])
            if residues not in alm:
                new_score = thisScore + 1
            else:
                new_score = thisScore
            fullScore = new_score + max(len(restA)-1,len(restB)-1)
            if fullScore <= maxScore:
                maxScore = fullScore
                queue += [[alm+[residues],new_score,restA[1:],restB[1:]]]

        # start adding gap
        if restA:
            residues = (restA[0],'-')
            if residues not in alm:
                new_score = thisScore + 1
            else:
                new_score = thisScore

            fullScore = new_score + max(len(restA)-1,len(restB))
            if fullScore <= maxScore:
                queue += [[alm+[residues],new_score,restA[1:],restB]]
        
        # add gap in a
        if restB:
            residues = ('-',restB[0])
            if residues not in alm:
                new_score = thisScore + 1
            else:
                new_score = thisScore

            fullScore = new_score + max(len(restA),len(restB)-1)
            if fullScore <= maxScore:
                queue += [[alm+[residues],new_score,restA,restB[1:]]]

    return out,maxScore


if __name__ == '__main__':
    test = structalign(argv[1],argv[2])
    for t in test[0]:
        print(' '.join(t[0]))
        print(' '.join(t[1]))
        print('-' * 2 * len(t[0]))
    print(test[1])
