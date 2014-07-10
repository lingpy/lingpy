from numpy import *
from scipy import stats
from lingpy import *

from itertools import product

lex = LexStat('samoyedic.csv')

langs = lex.language

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(langs)) + " languages: " + ", ".join(langs));

class Data:
    sigma0 = [] #alphabet of correspondences in the proto language model
    sigma1 = dict() #alphabet of segments in first observed language
    sigma2 = dict() #alphabet of segments in first observed language
    pairs = [] #cognate candidate pairs (lists of IPA tokens]
    
def build_data(lexstat, l1, l2, u, v):
    data = Data()
    #STEP 1: extract the candidate pairs from the LexStat object
    for concept in lex.concept[0:10]:
        entryIdxs = lex.get_list(concept=concept, flat=True)
        for iIdx in range(len(entryIdxs)):
            i = entryIdxs[iIdx]
            for jIdx in range(len(entryIdxs)):
                j = entryIdxs[jIdx]
                if lex[i][3] == l1 and lex[j][3] == l2:
                    word1 = lex[i][4]
                    word2 = lex[j][4]
                    word1.append("#")
                    word2.append("#")
                    data.pairs.append((word1, word2))
                    for segment in word1:
                        data.sigma1[segment] = True
                    for segment in word2:
                        data.sigma2[segment] = True
    print(str(len(data.pairs)) + " candidate pairs extracted for " + l1 + " and " + l2)
    print(l1 + " alphabet contains " + str(len(data.sigma1)) + " segments")
    print(l2 + " alphabet contains " + str(len(data.sigma2)) + " segments")
    #STEP 2: generate the alphabet of possible sound correspondences
    sequences1 = dict()
    sequences2 = dict()
    for seqLength in range(u,v+1):
        if seqLength == 0:
            sequences1[""] = True
            sequences2[""] = True
        elif seqLength == 1:
            for seq in data.sigma1.keys():
                sequences1[seq] = True
            for seq in data.sigma2.keys():
                sequences2[seq] = True 
        else:
            for (word1, word2) in data.pairs:
                for i in range(0,len(word1) - seqLength):
                    sequences1[".".join(word1[i:i+seqLength])] = True
                for i in range(0,len(word2) - seqLength):
                    sequences2[".".join(word2[i:i+seqLength])] = True
    data.sigma0 = list(product(sequences1, sequences2))
    #print(data.sigma0)
    print("(" + str(u) + "," + str(v) + ")-Model contains " + str(len(data.sigma0)) + " possible correspondences.")
    return data

class Hypothesis:
    f0 = dict() #generative model of the proto language, probabilities indexed by correspondences
    f1 = dict() #generative model of the first observed language, probabilities indexed by segments
    f2 = dict() #generative model of the second observed language, probabilities indexed by segments
    c = dict() #cognacy judgments (values between 0 and 1 for each pair)

def build_random_hypothesis(data):
    #TODO: sample dirichlet distribution without pseudo-observations?
    #for now: start with uniform distribution over f parameters, random cognacy judgments
    hypothesis = Hypothesis()
    hypothesis.f0 = dict([(symbol, 1.0 / size(data.sigma0)) for symbol in data.sigma0])
    hypothesis.f1 = dict([(symbol, 1.0 / size(data.sigma1)) for symbol in data.sigma1.keys()])
    hypothesis.f2 = dict([(symbol, 1.0 / size(data.sigma2)) for symbol in data.sigma2.keys()])
    hypothesis.c = dict([("".join(pair[0]) + "\t" + "".join(pair[1]), random.random()) for pair in data.pairs])
    return hypothesis
    
def eval_likelihood(data, hypothesis):
    #TODO: implement the expression for P(D|h) derived by Ellison (2007)
    logScore = 0.0
    for (word1, word2) in data.pairs:
        cognateScore = 0.0 #TODO: implement this
        word1Score = word_probability(word1, hypothesis.f1)
        word2Score = word_probability(word2, hypothesis.f2)
        nonCognateScore = word1Score * word2Score #non-cognates are generated independently
        cogProb = hypothesis.c["".join(word1) + "\t" + "".join(word2)]
        pairScore = cognateScore * cogProb + nonCognateScore * (1.0 - cogProb)
        logScore += log(pairScore) #product over all pairs
    return exp(logScore)

#use this to compute P1 and P2
def word_probability(word, alphabet_dist):
    score = 1.0
    for symbol in word:
        score *= alphabet_dist[symbol]
    return score

def gradient_ascent_step(data, hypothesis):
    #TODO: compute all derivatives of the likelihood function
    #TODO: implement the matrix and use it to generate the next hypothesis
    pass

#here comes the test program
data = build_data(lex, "Nenets", "Nganasan", 0,1)
hypothesis = build_random_hypothesis(data)
print(eval_likelihood(data, hypothesis))

#alternative that generates all sequences:
#         if seqLength == 0:
#             sequences1[""] = True
#             sequences2[""] = True
#         elif seqLength == 1:
#             for seq in data.sigma1.keys():
#                 sequences1[seq] = True
#             for seq in data.sigma2.keys():
#                 sequences2[seq] = True   
#         else:
#             for seq in product(data.sigma1.keys(), data.sigma1.keys(), repeat = seqLength - 1):
#                 sequences1[".".join(seq)] = True
#             for seq in product(data.sigma2.keys(), data.sigma2.keys(), repeat = seqLength - 1):
#                 sequences2[".".join(seq)] = True