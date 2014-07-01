from numpy import *
from lingpy import *

lex = LexStat('../data/eurasian/samoyedic.csv')

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(lex.language)) + " languages: " + ", ".join(lex.language));

lex.get_scorer()

print("Positive correspondences in scorer: ")
invertedMap = {v:k for k, v in lex.rscorer.chars2int.items()}
dim = len(lex.rscorer.matrix)
for i in range(dim):
    for j in range(dim):
        if (lex.rscorer.matrix[i][j] > 0.0):
            print(invertedMap[i] + "<->" + invertedMap[j] + ": " + str(lex.rscorer.matrix[i][j]))

print("\nCorrespondences across Dolgopolsky classes:")        
for i in range(dim):
    for j in range(dim):
        if lex.rscorer.matrix[i][j] > 0.0:
            class1 = invertedMap[i][0]
            class2 = invertedMap[j][0]
            if class1 != class2:
                print(invertedMap[i] + "<->" + invertedMap[j] + ": " + str(lex.rscorer.matrix[i][j]))
                
print("\nStability of classes at identical prosodic positions:")        
for i in range(dim):
    print(invertedMap[i] + ":" + str(lex.rscorer.matrix[i][i]))