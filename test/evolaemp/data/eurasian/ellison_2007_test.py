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
    for concept in lex.concept[0:100]:
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
            for seq1 in data.sigma1.keys():
                sequences1[seq1] = True
            for seq2 in data.sigma2.keys():
                sequences2[seq2] = True            
        else:
            for (word1, word2) in data.pairs:
                for i in range(0,len(word1) - seqLength):
                    sequences1[".".join(word1[i:i+seqLength])] = True
                for j in range(0,len(word2) - seqLength):
                    sequences2[".".join(word2[j:j+seqLength])] = True
    #data.sigma0 = list(product(sequences1, sequences2))
    pairs0 = dict()
    for seqLength1 in range(u,v+1):
        for seqLength2 in range(u,v+1):
            for (word1, word2) in data.pairs:
                for i in range(0,len(word1) - seqLength1 + 1):
                    for j in range(0,len(word2) - seqLength2 + 1):
                        pairs0[(".".join(word1[i:i+seqLength1]),".".join(word2[j:j+seqLength2]))] = True  
    data.sigma0 = pairs0.keys()
    #print(data.sigma0)
    print("(" + str(u) + "," + str(v) + ")-Model contains " + str(len(data.sigma0)) + " possible correspondences.")
    return data

#mapping variables to their indices in the hypothesis vector used during optimization
class HypothesisVariableIndex:
    f0 = dict() #generative model of the proto language, idexed by correspondences
    f1 = dict() #generative model of the first observed language, indexed by segments
    f2 = dict() #generative model of the second observed language, indexed by segments
    c = dict() #cognacy judgment variables
    #f0start, f1start, f2start, cstart for the additional indices

def build_random_hypothesis(data):
    #FIRST STEP: build the variable index
    hypothesis_index = HypothesisVariableIndex()
    varIndex = 0
    hypothesis_index.f0start = varIndex
    for symbol in data.sigma0:
        hypothesis_index.f0[symbol] = varIndex
        varIndex += 1
    hypothesis_index.f1start = varIndex
    for symbol in data.sigma1:
        hypothesis_index.f1[symbol] = varIndex
        varIndex += 1
    hypothesis_index.f2start = varIndex
    for symbol in data.sigma2:
        hypothesis_index.f2[symbol] = varIndex
        varIndex += 1
    hypothesis_index.cstart = varIndex   
    for pair in data.pairs:
        hypothesis_index.c["".join(pair[0]) + "\t" + "".join(pair[1])] = varIndex
        varIndex += 1
    #TODO: sample dirichlet distribution without pseudo-observations?
    #TODO: initialize f1 and f2 with observed segment distributions, c with string-based similarity?
    #for now: start with uniform distribution over f parameters, random cognacy judgments
    hypothesis_vector = zeros(varIndex)
    uniform_prob_f0 = 1.0 / len(data.sigma0)
    for index in hypothesis_index.f0.values():
        hypothesis_vector[index] = uniform_prob_f0
    uniform_prob_f1 = 1.0 / len(data.sigma1)
    for index in hypothesis_index.f1.values():
        hypothesis_vector[index] = uniform_prob_f1
    uniform_prob_f2 = 1.0 / len(data.sigma2)
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
            if pairScore == 0.0:
                #print("logScore = " + str(logScore))
                print("likelihood = " + str(0.0))
                return 0.0
            logScore += log(pairScore) #product over all pairs
        print("logLikelihood = " + str(logScore))
        #print("likelihood = " + str(exp(logScore)))
        return -logScore #exp(logScore) might have been to small!
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
    print("Most likely sound correspondences:")
    def key_func(key):
        #print("key_func(" + str(key) + " -> " + str(index.f0[key]) + ") = " + str(hypothesis[index.f0[key]]))
        return hypothesis[index.f0[key]]
    for pair in sorted(index.f0.keys(),key=key_func):
        print(str(pair) + "\t" + str(hypothesis[index.f0[pair]]))
    print("Cognate scores:")
    for pair in sorted(index.c.keys(), key= lambda x: hypothesis[index.c[x]]):
        print(pair + "\t" + str(hypothesis[index.c[pair]]))

#here comes the test program
data = build_data(lex, "Nenets", "Nganasan", 0,1)
for i in range(10):
    (index, hypothesis) = build_random_hypothesis(data)
    f = make_eval_likelihood(data, index)
    print(f(hypothesis))
equality_constraints = []
equality_constraints.append(lambda hypothesis : sum(hypothesis[0:index.f1start]) - 1)
equality_constraints.append(lambda hypothesis : sum(hypothesis[index.f1start:index.f2start]) - 1)
equality_constraints.append(lambda hypothesis : sum(hypothesis[index.f2start:index.cstart]) - 1)
result = optimize.fmin_slsqp(f, hypothesis, eqcons=equality_constraints, bounds=[(0,1) for i in range(size(hypothesis))])
print(result)
#result = optimize.fmin_cg(f, hypothesis)
print_hypothesis_summary(index, result)
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