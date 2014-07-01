from numpy import *
from lingpy import *

lex = LexStat('samoyedic.csv')

langs = lex.language

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(langs)) + " languages: " + ", ".join(langs));

lex.get_scorer(preprocessing=True)

for concept in lex.concept:
    dictset = dict()
    entryIdxs = lex.get_list(concept=concept, flat=True)
    for iIdx in range(len(entryIdxs)):
        i = entryIdxs[iIdx]
        cognateFound = False
        for jIdx in range(iIdx + 1, len(entryIdxs)):
            j = entryIdxs[jIdx]
            if lex[i][3] != lex[j][3]: #different languages
                if lex.align_pairs(i,j,return_distance=True,pprint=False) <= 0.5:
                    cognateFound = True
                    dictset[lex[i][3]] = lex[i][1]
                    dictset[lex[j][3]] = lex[j][1]
                    #print(concept + "\t " + lex[i][1] + " in " + lex[i][3] + ", " + lex[j][1] + " in " + lex[j][3])
        if not cognateFound and not lex[i][3] in dictset.keys():
            dictset[lex[i][3]] = "(" + lex[i][1] + ")"
    print(concept + "\t" + "\t".join([dictset[lang] for lang in lex.language]))
        
                    


