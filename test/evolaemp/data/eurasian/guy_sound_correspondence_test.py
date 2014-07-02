from numpy import *
from lingpy import *

from collections import defaultdict

import functools

lex = LexStat('samoyedic.csv')

langs = lex.language

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(langs)) + " languages: " + ", ".join(langs));

#lex.get_scorer(preprocessing=True)

def sum_row(matrix, rowIndex):
    return sum([value for value in matrix[rowIndex].values()])

def sum_column(matrix, columnIndex):
    return sum([matrix[rowIndex][columnIndex] for rowIndex in matrix.keys()])

def sum_matrix(matrix):
    return sum([sum_row(matrix, rowIndex) for rowIndex in matrix.keys()])

matrixA = defaultdict(lambda:defaultdict(lambda: 0))

for concept in lex.concept:
    entryIdxs = lex.get_list(concept=concept, flat=True)
    for iIdx in range(len(entryIdxs)):
        i = entryIdxs[iIdx]
        for jIdx in range(len(entryIdxs)):
            j = entryIdxs[jIdx]
            if lex[i][3] == 'Nenets' and lex[j][3] == 'Nganasan': #example Nenets vs. Nganasan
                smallerLength = min(len(lex[i][4]),len(lex[j][4]))
                for k in range(smallerLength):
                    matrixA[lex[i][4][k]][lex[j][4][k]] += 1
                    
langXKeys = matrixA.keys()
yKeyLists = [set(matrixA[xKey].keys()) for xKey in langXKeys]

langYKeys = functools.reduce(set.union,yKeyLists)

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixA[x][y]) for y in langYKeys]))

N = sum_matrix(matrixA)

print("N:" + str(N))
    
matrixB = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        matrixB[x][y] = (sum_row(matrixA,x)*sum_column(matrixA,y))/N
        
print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(int(round(matrixB[x][y],0))) for y in langYKeys]))
        
matrixE = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        matrixE[x][y] = matrixA[x][y] - int(round(matrixB[x][y],0))

print("\n")

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixE[x][y]) for y in langYKeys]))

scores = dict(((x,y),matrixE[x][y]) for x in langXKeys for y in langYKeys)
      
for pair in sorted(scores, key=scores.get):
    if scores[pair] > 19 or scores[pair] < -19:
        print(str(pair) + ": " + str(scores[pair]))