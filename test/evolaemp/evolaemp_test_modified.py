from numpy import *
from lingpy2 import *
import math

#READING IN THE ASJP DATA

f = open('data/asjp/sounds41.txt','r')
sounds = array([x.strip() for x  in f.readlines()])
f.close()

# reading in 'asjpMatrix.txt' into a numpy array
mfile = open("data/asjp/asjpMatrix.txt","r")
asjpRaw = mfile.readlines()
mfile.close()

asjpMatrix = array([x.strip().split('\t') for x in asjpRaw])


# restricting asjpMatrix to the languages in world_names.txt
f = open('data/asjp/world_names.txt','r')
rl = f.readlines()
f.close()

names = array([x.strip() for x in rl])

asjpMatrix = array([x for x in asjpMatrix if x[0] in names])

#INITIALZING LINGPY

# get the evolaemp schema
rc(schmema='evolaemp')

#TEST 1: LINGPY-BASED LANGUAGE DISTANCE MEASURE
print("\nTest 1: language distance measure based on pairwise alignment")
print("----------------------------------------------------------")

# Lingpy-based distance between row l1 and row l2 in mtx
def ldistLingpy(l1,l2,mtx=asjpMatrix):
    l1List = mtx[l1,4:]
    l2List = mtx[l2,4:]
    pairs = []
    wordOffsets = []
    for wordID in range(0,len(l1List)):
      if l1List[wordID] != '0' and l2List[wordID] != '0':
        wordOffsets.append(len(pairs))
        for entry1 in l1List[wordID].split('-'):
          for entry2 in l2List[wordID].split('-'):
            pairs.append((entry1,entry2))
    wordOffsets.append(len(pairs))
    #align all word pairs in parallel
    pair = Pairwise(pairs,merge_vowels=False)
    pair.align(distance=True,model=rcParams['asjp'])
    #collect the lowest distance values for each wordID (i.e. concept)
    distValues = [min(alignment[2] for alignment in pair.alignments[wordOffsets[i]:wordOffsets[i+1]]) for i in range(0,len(wordOffsets)-1)]
    #simply compute the average distance value
    averageDist = sum(distValues) / len(distValues)
    print("  ldistLingpy(" + asjpMatrix[l1][0] + "," + asjpMatrix[l2][0] + ") = " + str(averageDist))

ldistLingpy(1,2)
ldistLingpy(3,10)
ldistLingpy(5,6)
ldistLingpy(120,150)
ldistLingpy(235,1810)

#TEST 2: EXTRACTING PHONEME REPLACEMENT COUNTS USING MSA
print("\nTest 2: Extracting Phoneme Replacement Counts Using MSA")
print("----------------------------------------------------------")

lexdict = {}
lexdict[0] = ["ID", "concept", "ipa", "doculect"]
ID = 1

langs = range(1426,1460) #germanic languages

replacementOccurrenceTable = dict((asjpMatrix[taxon1,0],{}) for taxon1 in langs)
for taxon1 in langs:
  for taxon2 in langs:
    replacementOccurrenceTable[asjpMatrix[taxon1,0]][asjpMatrix[taxon2,0]] = {} 

#create a dictionary for cognate detection
for langID in langs:
  entries = asjpMatrix[langID,39] #entry for "mountain"
  for entry in entries.split('-'):
    lexdict[ID] = [langID, "mountain", entry, asjpMatrix[langID,0]]
    ID += 1

#cluster words into cognate sets
# XXX note that merge_vowels is now automatically set when choosing evolaemp as
# schema JML XXX
lexstat = LexStat(lexdict,model=rcParams['asjp'])
# this does not really work with only one entry, use the other sca-method
# instead
lexstat.cluster(method='sca',threshold=0.6)
etym_dict = lexstat.get_etymdict(ref='scaid', entry='', loans=False)

for cognateID in etym_dict.keys():
    entry_msq_file = open("cognate" + str(cognateID) + ".msq", 'w')
    entry_msq_file.write("ASJP database\n")
    entry_msq_file.write("Cognate " + str(cognateID) + " for Germanic languages\n")
    for IDList in etym_dict[cognateID]:
      if (IDList != 0):
        [langID, word, entry, langName] = lexdict[IDList[0]][:4]
        entry_msq_file.write(langName + "\t" + entry + "\n")
    entry_msq_file.close()
    print("Aligning cognate " + str(cognateID) + ":\n")
    multi = MSA("./cognate" + str(cognateID) + ".msq",merge_vowels=False)
    multi.prog_align(model=rcParams['asjp'],gop=-2,scale=0.7,factor=0.3)
    print(multi)
    #collect the sound replacements in this cognate set
    cognateSize = len(multi.seqs)
    if cognateSize > 1:
      mtx = multi.alm_matrix
      length = len(mtx[0])
      for i in range(0,cognateSize):
        for j in range(0,cognateSize):
          for k in range(0,length):
            #print(multi.taxa[i] + "\t->\t" + multi.taxa[j] + ":\t" + mtx[i][k] + "/" + mtx[j][k])
            occurrenceDict = replacementOccurrenceTable[multi.taxa[i]][multi.taxa[j]]
            #print(occurrenceDict)
            if mtx[i][k] not in occurrenceDict:
              occurrenceDict[mtx[i][k]] = {}
            if mtx[j][k] not in occurrenceDict[mtx[i][k]]:
              occurrenceDict[mtx[i][k]][mtx[j][k]] = 0.0
            occurrenceDict[mtx[i][k]][mtx[j][k]] += 1.0
            
#normalize the numbers in the occurrence dictionary
for taxon1 in replacementOccurrenceTable.keys():
  for taxon2 in replacementOccurrenceTable[taxon1].keys():
    for phon1 in replacementOccurrenceTable[taxon1][taxon2].keys():
      entries = replacementOccurrenceTable[taxon1][taxon2][phon1]
      entrySum = sum(list(entries.values()))
      for phon2 in entries.keys():
        print("  " + taxon1 + "->" + taxon2 + ", " + phon1 + "->" + phon2 + ": " + str(entries[phon2]) + "/" + str(entrySum))
        entries[phon2] = entries[phon2] / entrySum
