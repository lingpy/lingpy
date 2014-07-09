from numpy import *
from scipy import stats
from lingpy import *

from collections import defaultdict

import functools
import math

scoreComputation = "z"  # diff, z, chi2, G2

lex = LexStat('samoyedic.csv')

langs = lex.language

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(langs)) + " languages: " + ", ".join(langs));

# lex.get_scorer(preprocessing=True)

def sum_row(matrix, rowIndex):
    return sum([value for value in matrix[rowIndex].values()])

def sum_column(matrix, columnIndex):
    return sum([matrix[rowIndex][columnIndex] for rowIndex in matrix.keys()])

def sum_matrix(matrix):
    return sum([sum_row(matrix, rowIndex) for rowIndex in matrix.keys()])

def sum_reduced_matrix(matrix, rowIndex, columnIndex):
    return sum_matrix(matrix) - sum_row(matrix, rowIndex) - sum_column(matrix, columnIndex) + matrix[rowIndex][columnIndex]

def row_variance(matrix, rowIndex):
    expected = sum_row(matrix, rowIndex) / len(matrix[rowIndex])
    return sum([(value - expected) ** 2 for value in matrix[rowIndex].values()]) / len(matrix[rowIndex])

def contract_alignments(word1, word2):
    for i in range(len(word1)-1):
        if word1[i] == "-" == word2[i+1]:
            del word1[i]
            del word2[i+1]
        elif word1[i+1] == "-" == word2[i]:
            del word1[i+1]
            del word2[i]   

def print_alignment_details(word1, word2, weight_matrix):
    print(" " + "\t".join(word1))
    print(" " + "\t".join(word2))  
    print(" " + "\t".join(map(str,map(int,[weight_matrix[word1[i]][word2[i]] for i in range(len(word1))]))))   

matrixA = defaultdict(lambda:defaultdict(lambda: 0))

# FIRST PHASE: COUNT CO-OCCURRENCES
for concept in lex.concept:
    entryIdxs = lex.get_list(concept=concept, flat=True)
    for iIdx in range(len(entryIdxs)):
        i = entryIdxs[iIdx]
        for jIdx in range(len(entryIdxs)):
            j = entryIdxs[jIdx]
            if lex[i][3] == 'Nenets' and lex[j][3] == 'Nganasan':  # example Nenets vs. Nganasan
                for k in range(len(lex[i][4])):
                    for l in range(len(lex[j][4])):
                        matrixA[lex[i][4][k]][lex[j][4][l]] += 1
                    
langXKeys = matrixA.keys()
yKeyLists = [set(matrixA[xKey].keys()) for xKey in langXKeys]

langYKeys = functools.reduce(set.union, yKeyLists)

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixA[x][y]) for y in langYKeys]))

N = sum_matrix(matrixA)

print("N:" + str(N))
    
matrixB = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        matrixB[x][y] = (sum_row(matrixA, x) * sum_column(matrixA, y)) / N
        
print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixB[x][y]) for y in langYKeys]))
        
matrixE = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        if matrixA[x][y] <= 20:
            matrixE[x][y] = 1 - stats.binom.sf(matrixA[x][y] - 1, matrixB[x][y], sum_column(matrixA, y) / N)
        else:    
            if scoreComputation == "diff":
                matrixE[x][y] = matrixA[x][y] - matrixB[x][y]  # simple difference: Observed - Expected
            elif scoreComputation == "z":
                # matrixE[x][y] = (matrixA[x][y] - matrixB[x][y])/math.sqrt(row_variance(matrixB, x)) #z score (??): (Observed - Expected)/standard deviation of expected count
                # pooled 2-sample z-score test
                zScoreAbove = matrixA[x][y] / sum_row(matrixA, x) - (sum_column(matrixA, y) - matrixA[x][y]) / (N - sum_row(matrixA, x))
                zScoreBelow = sum_column(matrixA, y) / N
                zScoreBelow = zScoreBelow * (1 - zScoreBelow)
                zScoreBelow /= (N * (sum_row(matrixA, x) / N) * (1 - sum_row(matrixA, x) / N))
                zScoreBelow = math.sqrt(zScoreBelow)
                matrixE[x][y] = zScoreAbove / zScoreBelow
            elif scoreComputation == "chi2":
                # matrixE[x][y] = (matrixA[x][y] - matrixB[x][y])**2/matrixB[x][y] #chi square: (Observed - Expected)^2/Expected
                chi2Sum = (matrixA[x][y] - matrixB[x][y]) ** 2 / matrixB[x][y]
                chi2Sum += (sum_reduced_matrix(matrixA, x, y) - sum_reduced_matrix(matrixB, x, y)) ** 2 / sum_reduced_matrix(matrixB, x, y)
                chi2Sum += ((sum_row(matrixA, x) - matrixA[x][y]) - (sum_row(matrixB, x) - matrixB[x][y])) ** 2 / (sum_row(matrixB, x) - matrixB[x][y])
                chi2Sum += ((sum_column(matrixA, y) - matrixA[x][y]) - (sum_column(matrixB, y) - matrixB[x][y])) ** 2 / (sum_column(matrixB, y) - matrixB[x][y])
                matrixE[x][y] = chi2Sum
            elif scoreComputation == "G2":
                if matrixA[x][y] == 0:
                    matrixE[x][y] = 0
                else:
                    matrixE[x][y] = matrixA[x][y] * math.log((matrixA[x][y] / matrixB[x][y]))  # G square: Observed * log(Observed/Expected)

print("\n")

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixE[x][y]) for y in langYKeys]))

# OPTIMIZATIONS USING THRESHOLD VALUES AS STATED IN THE PAPER
for x in langXKeys:
    for y in langYKeys:
        if matrixE[x][y] > 4.0:
            matrixE[x][y] = 4.0
        elif matrixE[x][y] > 0.0 and matrixE[x][y] < 1.28:
            matrixE[x][y] = 0.0
            
scores = dict(((x, y), matrixE[x][y]) for x in langXKeys for y in langYKeys)
      
for pair in sorted(scores, key=scores.get):
    if scores[pair] != 0.0:
        print(str(pair) + ": " + str(scores[pair]))

# COMPUTE THE ALIGNMENTS
cognacies = dict()
alignments = []
for concept in lex.concept:
    entryIdxs = lex.get_list(concept=concept, flat=True)
    for iIdx in range(len(entryIdxs)):
        i = entryIdxs[iIdx]
        for jIdx in range(len(entryIdxs)):             
           j = entryIdxs[jIdx]
           if lex[i][3] == 'Nenets' and lex[j][3] == 'Nganasan':
               word1 = lex[i][4]
               word2 = lex[j][4]
               #print("(" + "".join(word1) + "," + "".join(word2) + ")")
               # fill potential table
               pot = empty((len(word1), len(word2)))
               for k in reversed(range(len(word1))):
                   for l in reversed(range(len(word2))):
                       pot[k][l] = matrixE[word1[k]][word2[l]] 
                       if k < len(word1) - 1 and l < len(word2) - 1:
                           pot[k][l] += amax(pot[k + 1:, l + 1:])
               # compute the alignment
               word1Aligned = []
               word2Aligned = []
               pos1 = 0
               pos2 = 0
               while pos1 < len(word1) and pos2 < len(word2):
                   #print(pot[pos1:, pos2:])
                   #print("argmax= " + str(argmax(pot[pos1:, pos2:])))
                   (posDiff1, posDiff2) = unravel_index(argmax(pot[pos1:, pos2:]), (pot.shape[0] - pos1, pot.shape[1] - pos2))
                   for q in range(min(posDiff1,posDiff2)):
                       word1Aligned.append(word1[pos1 + q])
                       word2Aligned.append(word2[pos2 + q])
                   for q in range(max(posDiff1-posDiff2,0)):
                       word1Aligned.append(word1[pos1 + min(posDiff1,posDiff2) + q])  
                       word2Aligned.append('-')   
                   for q in range(max(posDiff2-posDiff1,0)):
                       word1Aligned.append('-')  
                       word2Aligned.append(word2[pos2 + min(posDiff1,posDiff2) + q])    
                   word1Aligned.append(word1[pos1+posDiff1])
                   word2Aligned.append(word2[pos2+posDiff2]) 
                   pos1 = pos1 + posDiff1 + 1
                   pos2 = pos2 + posDiff2 + 1
                   #print("(" + str(posDiff1) + "," + str(posDiff2) + ") -> (" + str(pos1) + "," + str(pos2) + ")")
               for q in range(max(len(word1) - pos1,0)):
                   word1Aligned.append(word1[pos1 + q])  
                   word2Aligned.append('-')  
               for q in range(max(len(word2) - pos2,0)):
                   word1Aligned.append('-')  
                   word2Aligned.append(word2[pos2 + q])   
               contract_alignments(word1Aligned,word2Aligned)
               alignments.append((word1Aligned,word2Aligned))
               print("(" + "".join(word1) + "," + "".join(word2) + "): " + str(amax(pot) / max(len(word1), len(word2))))
               print_alignment_details(word1Aligned, word2Aligned, matrixE)
               cognacies["(" + "".join(word1) + "," + "".join(word2) + ")"] = str(amax(pot) / max(len(word1), len(word2)))

for pair in sorted(cognacies, key=cognacies.get):
        # pass
        print(str(pair) + ":" + cognacies[pair])

#RESCORING BASED ON ALIGNMENTS  
matrixA = defaultdict(lambda:defaultdict(lambda: 0))
    
for (w1, w2) in alignments:
    for k in range(len(w1)):
        matrixA[w1[k]][w2[k]] += 1


#REPEAT THE ABOVE PROCESS (TODO: MODULARIZE)
                    
langXKeys = matrixA.keys()
yKeyLists = [set(matrixA[xKey].keys()) for xKey in langXKeys]

langYKeys = functools.reduce(set.union, yKeyLists)

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixA[x][y]) for y in langYKeys]))

N = sum_matrix(matrixA)

print("N:" + str(N))
    
matrixB = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        matrixB[x][y] = (sum_row(matrixA, x) * sum_column(matrixA, y)) / N
        
print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixB[x][y]) for y in langYKeys]))
        
matrixE = defaultdict(lambda:defaultdict(lambda: 0.0))
for x in langXKeys:
    for y in langYKeys:
        if matrixA[x][y] <= 20:
            matrixE[x][y] = 1 - stats.binom.sf(matrixA[x][y] - 1, matrixB[x][y], sum_column(matrixA, y) / N)
        else:    
            if scoreComputation == "diff":
                matrixE[x][y] = matrixA[x][y] - matrixB[x][y]  # simple difference: Observed - Expected
            elif scoreComputation == "z":
                # matrixE[x][y] = (matrixA[x][y] - matrixB[x][y])/math.sqrt(row_variance(matrixB, x)) #z score (??): (Observed - Expected)/standard deviation of expected count
                # pooled 2-sample z-score test
                zScoreAbove = matrixA[x][y] / sum_row(matrixA, x) - (sum_column(matrixA, y) - matrixA[x][y]) / (N - sum_row(matrixA, x))
                zScoreBelow = sum_column(matrixA, y) / N
                zScoreBelow = zScoreBelow * (1 - zScoreBelow)
                zScoreBelow /= (N * (sum_row(matrixA, x) / N) * (1 - sum_row(matrixA, x) / N))
                zScoreBelow = math.sqrt(zScoreBelow)
                matrixE[x][y] = zScoreAbove / zScoreBelow
            elif scoreComputation == "chi2":
                # matrixE[x][y] = (matrixA[x][y] - matrixB[x][y])**2/matrixB[x][y] #chi square: (Observed - Expected)^2/Expected
                chi2Sum = (matrixA[x][y] - matrixB[x][y]) ** 2 / matrixB[x][y]
                chi2Sum += (sum_reduced_matrix(matrixA, x, y) - sum_reduced_matrix(matrixB, x, y)) ** 2 / sum_reduced_matrix(matrixB, x, y)
                chi2Sum += ((sum_row(matrixA, x) - matrixA[x][y]) - (sum_row(matrixB, x) - matrixB[x][y])) ** 2 / (sum_row(matrixB, x) - matrixB[x][y])
                chi2Sum += ((sum_column(matrixA, y) - matrixA[x][y]) - (sum_column(matrixB, y) - matrixB[x][y])) ** 2 / (sum_column(matrixB, y) - matrixB[x][y])
                matrixE[x][y] = chi2Sum
            elif scoreComputation == "G2":
                if matrixA[x][y] == 0:
                    matrixE[x][y] = 0
                else:
                    matrixE[x][y] = matrixA[x][y] * math.log((matrixA[x][y] / matrixB[x][y]))  # G square: Observed * log(Observed/Expected)

print("\n")

print ("\t" + "\t".join(langYKeys))

for x in langXKeys:
    print(x + ":\t" + "\t".join([str(matrixE[x][y]) for y in langYKeys]))

# OPTIMIZATIONS USING THRESHOLD VALUES AS STATED IN THE PAPER
for x in langXKeys:
    for y in langYKeys:
        if matrixE[x][y] > 4.0:
            matrixE[x][y] = 4.0
        elif matrixE[x][y] > 0.0 and matrixE[x][y] < 1.28:
            matrixE[x][y] = 0.0
            
scores = dict(((x, y), matrixE[x][y]) for x in langXKeys for y in langYKeys)
      
for pair in sorted(scores, key=scores.get):
    if scores[pair] != 0.0:
        print(str(pair) + ": " + str(scores[pair]))

# COMPUTE THE ALIGNMENTS
cognacies = dict()
alignments = []
for concept in lex.concept:
    entryIdxs = lex.get_list(concept=concept, flat=True)
    for iIdx in range(len(entryIdxs)):
        i = entryIdxs[iIdx]
        for jIdx in range(len(entryIdxs)):             
           j = entryIdxs[jIdx]
           if lex[i][3] == 'Nenets' and lex[j][3] == 'Nganasan':
               word1 = lex[i][4]
               word2 = lex[j][4]
               #print("(" + "".join(word1) + "," + "".join(word2) + ")")
               # fill potential table
               pot = empty((len(word1), len(word2)))
               for k in reversed(range(len(word1))):
                   for l in reversed(range(len(word2))):
                       pot[k][l] = matrixE[word1[k]][word2[l]] 
                       if k < len(word1) - 1 and l < len(word2) - 1:
                           pot[k][l] += amax(pot[k + 1:, l + 1:])
               # compute the alignment
               word1Aligned = []
               word2Aligned = []
               pos1 = 0
               pos2 = 0
               while pos1 < len(word1) and pos2 < len(word2):
                   #print(pot[pos1:, pos2:])
                   #print("argmax= " + str(argmax(pot[pos1:, pos2:])))
                   (posDiff1, posDiff2) = unravel_index(argmax(pot[pos1:, pos2:]), (pot.shape[0] - pos1, pot.shape[1] - pos2))
                   for q in range(min(posDiff1,posDiff2)):
                       word1Aligned.append(word1[pos1 + q])
                       word2Aligned.append(word2[pos2 + q])
                   for q in range(max(posDiff1-posDiff2,0)):
                       word1Aligned.append(word1[pos1 + min(posDiff1,posDiff2) + q])  
                       word2Aligned.append('-')   
                   for q in range(max(posDiff2-posDiff1,0)):
                       word1Aligned.append('-')  
                       word2Aligned.append(word2[pos2 + min(posDiff1,posDiff2) + q])    
                   word1Aligned.append(word1[pos1+posDiff1])
                   word2Aligned.append(word2[pos2+posDiff2]) 
                   pos1 = pos1 + posDiff1 + 1
                   pos2 = pos2 + posDiff2 + 1
                   #print("(" + str(posDiff1) + "," + str(posDiff2) + ") -> (" + str(pos1) + "," + str(pos2) + ")")
               for q in range(max(len(word1) - pos1,0)):
                   word1Aligned.append(word1[pos1 + q])  
                   word2Aligned.append('-')  
               for q in range(max(len(word2) - pos2,0)):
                   word1Aligned.append('-')  
                   word2Aligned.append(word2[pos2 + q])   
               contract_alignments(word1Aligned,word2Aligned)
               alignments.append((word1Aligned,word2Aligned))
               print("(" + "".join(word1) + "," + "".join(word2) + "): " + str(amax(pot) / max(len(word1), len(word2))))
               print_alignment_details(word1Aligned, word2Aligned, matrixE)
               cognacies["(" + "".join(word1) + "," + "".join(word2) + ")"] = str(amax(pot) / max(len(word1), len(word2)))

for pair in sorted(cognacies, key=cognacies.get):
        # pass
        print(str(pair) + ":" + cognacies[pair])