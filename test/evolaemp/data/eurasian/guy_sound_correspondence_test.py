from numpy import *
from lingpy import *

from collections import defaultdict

import functools
import math

scoreComputation = "z" #diff, z, chi2, G2

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

def row_variance(matrix, rowIndex):
    expected = sum_row(matrix, rowIndex) / len(matrix[rowIndex])
    return sum([(value - expected)**2 for value in matrix[rowIndex].values()]) / len(matrix[rowIndex])

matrixA = defaultdict(lambda:defaultdict(lambda: 0))

#FIRST PHASE: COUNT CO-OCCURRENCES
for concept in lex.concept:
    entryIdxs = lex.get_list(concept=concept, flat=True)
    for iIdx in range(len(entryIdxs)):
        i = entryIdxs[iIdx]
        for jIdx in range(len(entryIdxs)):
            j = entryIdxs[jIdx]
            if lex[i][3] == 'Nenets' and lex[j][3] == 'Nganasan': #example Nenets vs. Nganasan
                for k in range(len(lex[i][4])):
                    for l in range(len(lex[j][4])):
                        matrixA[lex[i][4][k]][lex[j][4][l]] += 1
                    
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
    print(x + ":\t" + "\t".join([str(matrixB[x][y]) for y in langYKeys]))
        
matrixE = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        if scoreComputation == "diff":
            matrixE[x][y] = matrixA[x][y] - matrixB[x][y] #simple difference: Observed - Expected
        elif scoreComputation == "z":
            matrixE[x][y] = (matrixA[x][y] - matrixB[x][y])/math.sqrt(row_variance(matrixB, x)) #z score (??): (Observed - Expected)/standard deviation of expected count
        elif scoreComputation == "chi2":
            matrixE[x][y] = (matrixA[x][y] - matrixB[x][y])**2/matrixB[x][y] #chi square: (Observed - Expected)^2/Expected
        elif scoreComputation == "G2":
            if matrixA[x][y] == 0:
                matrixE[x][y] = 0
            else:
                matrixE[x][y] = matrixA[x][y] * math.log((matrixA[x][y]/matrixB[x][y])) #G square: Observed * log(Observed/Expected)

print("\n")

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixE[x][y]) for y in langYKeys]))

scores = dict(((x,y),matrixE[x][y]) for x in langXKeys for y in langYKeys)
      
for pair in sorted(scores, key=scores.get):
    if scores[pair] > 1.28 or scores[pair] < -0:
        print(str(pair) + ": " + str(scores[pair]))

# RESCORING BASED ON ALIGNMENTS      
# for concept in lex.concept:
#     entryIdxs = lex.get_list(concept=concept, flat=True)
#     for iIdx in range(len(entryIdxs)):
#         i = entryIdxs[iIdx]
#         for jIdx in range(len(entryIdxs)):
#             j = entryIdxs[jIdx]
#             if lex[i][3] == 'Nenets' and lex[j][3] == 'Nganasan': #example Nenets vs. Nganasan
#                 smallerLength = min(len(lex[i][4]),len(lex[j][4]))
#                 for k in range(smallerLength):
#                     matrixA[lex[i][4][k]][lex[j][4][k]] += 1