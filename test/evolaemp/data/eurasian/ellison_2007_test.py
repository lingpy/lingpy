from numpy import *
from scipy import stats
from scipy import optimize
from lingpy import *
import copy

from itertools import product

lex = LexStat('samoyedic.csv')

langs = lex.language

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(langs)) + " languages: " + ", ".join(langs));

verbose = False
spill_guts = False

class Data:
    sigma0 = [] #alphabet of correspondences in the proto language model
    sigma1 = dict() #alphabet of segments in first observed language
    sigma2 = dict() #alphabet of segments in first observed language
    pairs = [] #cognate candidate pairs (lists of IPA tokens]
    #added information: u and v (upper and lower boundaries for correspondence length
    
def build_data(lexstat, l1, l2, u, v):
    data = Data()
    data.u = u
    data.v = v
    #STEP 1: extract the candidate pairs from the LexStat object
    for concept in lex.concept[0:2]:
        entryIdxs = lex.get_list(concept=concept, flat=True)
        for iIdx in range(len(entryIdxs)):
            i = entryIdxs[iIdx]
            for jIdx in range(len(entryIdxs)):
                j = entryIdxs[jIdx]
                if lex[i][3] == l1 and lex[j][3] == l2:
                    word1 = copy.copy(lex[i][4])
                    word2 = copy.copy(lex[j][4])
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

#mapping variables to their indices in the hypothesis vector used during optimization
class HypothesisVariableIndex:
    f0 = dict() #generative model of the proto language, ndexed by correspondences
    f1 = dict() #generative model of the first observed language, indexed by segments
    f2 = dict() #generative model of the second observed language, indexed by segments
    c = dict() #cognacy judgment variables

def build_random_hypothesis(data):
    #FIRST STEP: build the variable index
    hypothesis_index = HypothesisVariableIndex()
    varIndex = 0
    for symbol in data.sigma0:
        hypothesis_index.f0[symbol] = varIndex
        varIndex += 1
    for symbol in data.sigma1:
        hypothesis_index.f1[symbol] = varIndex
        varIndex += 1
    for symbol in data.sigma2:
        hypothesis_index.f2[symbol] = varIndex
        varIndex += 1
    for pair in data.pairs:
        hypothesis_index.c["".join(pair[0]) + "\t" + "".join(pair[1])] = varIndex
        varIndex += 1
    #TODO: sample dirichlet distribution without pseudo-observations?
    #TODO: initialize f1 and f2 with observed segment distributions, c with string-based similarity?
    #for now: start with uniform distribution over f parameters, random cognacy judgments
    hypothesis_vector = zeros(varIndex)
    uniform_prob_f0 = 1.0 / size(data.sigma0)
    for index in hypothesis_index.f0.values():
        hypothesis_vector[index] = uniform_prob_f0
    uniform_prob_f1 = 1.0 / size(data.sigma1)
    for index in hypothesis_index.f1.values():
        hypothesis_vector[index] = uniform_prob_f1
    uniform_prob_f2 = 1.0 / size(data.sigma2)
    for index in hypothesis_index.f2.values():
        hypothesis_vector[index] = uniform_prob_f2
    for index in hypothesis_index.c.values():
        hypothesis_vector[index] = random.random()
    return (hypothesis_index, hypothesis_vector)

def make_eval_likelihood(data, index):     
    def eval_likelihood(hypothesis_vector):
        logScore = 0.0
        for (word1, word2) in data.pairs:          
            cognateScore = word_pair_probability(data.u, data.v, word1, word2, 0, 0, index.f0, hypothesis_vector, dict())
            if verbose:
                print("word_pair_probability(" + "".join(word1) + "," + "".join(word2) + ") = " + str(cognateScore))
            
            word1Score = word_probability(word1, index.f1, hypothesis_vector)
            word2Score = word_probability(word2, index.f2, hypothesis_vector)
            nonCognateScore = word1Score * word2Score #non-cognates are generated independently
            
            cogProb = hypothesis_vector[index.c["".join(word1) + "\t" + "".join(word2)]]
            pairScore = cognateScore * cogProb + nonCognateScore * (1.0 - cogProb)
            logScore += log(pairScore) #product over all pairs
        print("logScore = " + str(logScore))
        return -exp(logScore)
    return eval_likelihood

#use this to compute P1 and P2
def word_probability(word, alphabet_index, hypothesis):
    score = 1.0
    for symbol in word:
        score *= hypothesis[alphabet_index[symbol]]
    return score

def get_prefixes(u,v,word,startIndex):
    prefixes = []
    for length in range(u,min(v+1,len(word)-startIndex)):
        prefixes.append(word[startIndex:startIndex+length])
    return prefixes
    
#u and v: minimal and maximal prefix length
#word1 and word2: lists of tokens representing a cognate candidate pair
#startIndex1, startIndex: cursor positions (used in recursion to avoid creating list copies)
#correspondence_index: map correspondence pairs to the indexes of their probabilities in the hypothesis
#hypothesis: the hypothesis vector under which the word pair probability is evaluated
def word_pair_probability(u, v, word1, word2, startIndex1, startIndex2, correspondence_index, hypothesis, intermediate_results):
    if spill_guts:
        print("word_pair_probability(" + "".join(word1[startIndex1:]) + "," + "".join(word2[startIndex2:]) + ") = ")
    totalProbability = intermediate_results.get(("".join(word1[startIndex1:]),"".join(word2[startIndex2:])),None)
    if totalProbability != None:
        if spill_guts:
            print("(retrieved) " + str(totalProbability))
        return totalProbability
    if startIndex1 == len(word1)-1 and startIndex2 == len(word2)-1:
        if spill_guts:
            print("  hypothesis[" + str(correspondence_index[("#","#")]) + " = (#,#)] = " + str(hypothesis[correspondence_index[("#","#")]]))
        return hypothesis[correspondence_index[("#","#")]]
    else:
        totalProbability = 0
        for prefix1 in get_prefixes(u,v,word1,startIndex1):
            for prefix2 in get_prefixes(u,v,word2,startIndex2):
                if len(prefix1) == 0 and len(prefix2) == 0:
                    continue #otherwise infinite recursion
                else:
                    branchProbability = hypothesis[correspondence_index[("".join(prefix1),"".join(prefix2))]]
                    if spill_guts:        
                        print("  hypothesis[" + str(correspondence_index[("".join(prefix1),"".join(prefix2))]) + " = (" + "".join(prefix1) + "," + "".join(prefix2) + ")] = " + str(branchProbability))
                    branchProbability *= word_pair_probability(u, v, word1, word2, startIndex1 + len(prefix1), startIndex2 + len(prefix2), correspondence_index, hypothesis, intermediate_results)
                    totalProbability += branchProbability
        intermediate_results[("".join(word1[startIndex1:]),"".join(word2[startIndex2:]))] = totalProbability
        if spill_guts:
            print(str(totalProbability))
        return totalProbability
        
def print_hypothesis_summary(index, hypothesis):
    print("Cognate scores:")
    for pair in index.c.keys():
        print(pair + "\t" + str(hypothesis[index.c[pair]]))

#here comes the test program
data = build_data(lex, "Nenets", "Nganasan", 0,1)
for i in range(10):
    (index, hypothesis) = build_random_hypothesis(data)
    f = make_eval_likelihood(data, index)
    print(f(hypothesis))
print([(0,1) for i in range(size(hypothesis))])
result = optimize.fmin_l_bfgs_b(f, hypothesis, approx_grad=True, bounds=[(0,1) for i in range(size(hypothesis))])
print(result[0])
#result = optimize.fmin_cg(f, hypothesis)
print_hypothesis_summary(index, result[0])
#print(f(result))

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