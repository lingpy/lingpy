# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-12 11:56
# modified : 2013-03-12 11:56
"""
LexStat algorithm for automatic cognate detection.
"""

__author__="Johann-Mattis List"
__date__="2013-03-12"

# builtin
import random

# thirdparty modules
from ..thirdparty import cogent as cg

# lingpy-modules
from ..data import *
from ..sequence.sound_classes import *
from ..basic import Wordlist
try:
    from ..algorithm.cython import calign
except:
    from ..algorithm.cython import _calign as calign

class tdict(object):
    """

    """
    def __init__(self,chars,matrix):
        
        self.chars2int = dict([(char,i) for char,i in
            zip(chars,range(len(chars)))])
        self.matrix = matrix

    def __getitem__(self,x):
        
        try:
            return self.matrix[self.chars2int[x[0]]][self.chars2int[x[1]]]
        except:
            return -10.0

    def __repr__(self):
        return str(list(self.chars2int.keys()))

    def __str__(self):
        return str(list(self.chars2int.items()))

class MCBasic(object):

    def __init__(
            self,
            seqs
            ):
        """
        Basic object for Markov chains.
        """

        self.seqs = seqs
        
        # create distribution
        self.dist = {}
        for seq in self.seqs:

            for s1,s2 in zip(['#']+seq,seq+['$']):
                try:
                    self.dist[s1] += [s2]
                except:
                    self.dist[s1] = [s2]

        # create probabilities
        #self.probs = {}
        #for key,value in self.dist.items():
        #    # iterate over all items in the list
        #    for v in set(value):
        #        try:
        #            self.probs[key][v] = np.log10(self.dist[key].count(v)/len(self.dist[key]))
        #        except:
        #            self.probs[key] = {}
        #            self.probs[key][v] = np.log10(self.dist[key].count(v) / len(self.dist[key]))

    def walk(self):
        """
        Create random sequence from the distribution.
        """
        # out sequence
        out = []

        # get the start sequence 
        startS = random.choice(self.dist['#'])

        # add start to out
        out += [startS]

        # start looping
        while True:

            # get nextS
            nextS = random.choice(self.dist[out[-1]])

            # check for terminal symbol
            if nextS == '$':
                break
            
            out += [nextS]

        return out

    #def evaluate(
    #        self,
    #        seq
    #        ):
    #    """
    #    Evaluate the probability of a sequence to be created.
    #    """
    #    
    #    # probability set to zero
    #    prob = 0.0
    #    
    #    before = '#'
    #    for s in seq:
    #        try:
    #            prob += self.probs[before][s]
    #            before = s
    #        except:
    #            prob = 0.0
    #            break
    #    
    #    return prob
class MCPhon(MCBasic):
    
    def __init__(
            self,
            words,
            tokens=False,
            prostrings=[]
            ):
        """
        Markov Chains for phonetic sequences.
        """
        
        self.words = words
        self.tokens = []
        self.bigrams = []

        # start filling the dictionary
        for i,w in enumerate(words):
            
            # check for tokenized string
            if not tokens:
                tokens = ipa2tokens(w)
            else:
                tokens = w[:]
            self.tokens += [tokens]

            # create prosodic string
            if prostrings:
                p = prostrings[i]
            else:
                p = prosodic_string(tokens2class(tokens,model=art))

            # zip the stuff
            bigrams = list(zip(p,tokens))
            
            # start appending the stuff
            self.bigrams += [bigrams]

            # init the mother object
            MCBasic.__init__(self,self.bigrams)

    def get_string(
            self,
            new=True,
            tokens=False
            ):
        """
        Function generates a string based on the input data.
        """
        
        # create the first string
        out = self.walk()

        while new:
            if out in self.bigrams:
                out = self.walk()
            else:
                break

        if tokens:
            return out
        else:
            return ' '.join([i[1] for i in  out])

class LexStat(Wordlist):
    """

    """
    
    def __init__(
            self,
            filename,
            **keywords
            ):

        defaults = {
                "model" : sca
                }
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]
        
        # store the model
        self.model = keywords['model']

        # initialize the wordlist
        Wordlist.__init__(self,filename)
        
        # check for basic input data
        # tokens
        if not "tokens" in self.header:
            self.add_entries("tokens","ipa",lambda x:ipa2tokens(x))
        
        # sonority profiles
        if not "sonars" in self.header:
            self.add_entries(
                    "sonars",
                    "tokens",
                    lambda x:[int(i) for i in tokens2class(x,art)]
                    )

        # get prosodic strings
        if not "prostrings" in self.header:
            self.add_entries(
                    "prostrings",
                    "sonars",
                    lambda x:prosodic_string(x)
                    )
        
        # get sound class strings
        if not "classes" in self.header:
            self.add_entries(
                    "classes",
                    "tokens",
                    lambda x:''.join(tokens2class(x,keywords["model"]))
                    )
        
        # create IDs for the languages
        if not "langid" in self.header:
            transform = dict(zip(self.taxa,[str(i+1) for i in range(self.width)]))
            self.add_entries(
                    "langid",
                    "taxa",
                    lambda x:transform[x]
                    )
        # get the numbers for all strings
        if not "numbers" in self.header:
            self.add_entries(
                    "numbers",
                    "langid,classes,prostrings",
                    lambda x,y: ["{0}.{1[0]}.{1[1]}".format(x[y[0]],a) for a in zip(x[y[1]],x[y[2]])]    
                    )

        # check for weights
        if not "weights" in self.header:
            self.add_entries(
                    "weights",
                    "prostrings",
                    lambda x:prosodic_weights(x)
                    )

        # check for duplicates
        if not "duplicates" in self.header:
            duplicates = {}
            for taxon in self.taxa:
                words = []
                for idx in self.get_list(
                        col=taxon,
                        flat=True
                        ):
                    # get the words
                    word = self[idx,'words']
                    if word in words:
                        duplicates[idx] = 1
                    else:
                        duplicates[idx] = 0
                        words += [word]
            self.add_entries(
                    "duplicates",
                    duplicates,
                    lambda x:x
                    )

        # create an index 
        if not hasattr(self,'freqs'):
            self.chars = []
            self.freqs = {}
            for taxon in self.taxa:
                self.freqs[taxon] = {}
                words = self.get_list(
                        col=taxon,
                        entry='numbers',
                        flat=True
                        )
                for word in words:
                    for char in word:
                        try:
                            self.freqs[taxon][char] += 1
                        except:
                            self.freqs[taxon][char] = 1
                        self.chars.append(char)
            self.chars = list(set(self.chars))

        # create a scoring dictionary
        if not hasattr(self,"scorer"):
            matrix = [[0.0 for i in range(len(self.chars))] for j in range(len(self.chars))]
            for i,charA in enumerate(self.chars):
                for j,charB in enumerate(self.chars):
                    if i < j:
                        
                        # add dictionary scores to the scoredict
                        score = keywords["model"](
                                charA[charA.index('.')+1][0],
                                charB[charB.index('.')+1][0]
                                )
                        matrix[i][j] = score
                        matrix[j][i] = score
                    elif i == j:
                        # add dictionary scores to the scoredict
                        score = keywords["model"](
                                charA[charA.index('.')+1][0],
                                charB[charB.index('.')+1][0]
                                )
                        matrix[i][j] = score
        
            self.scorer = tdict(self.chars,matrix)

        # make the language pairs
        if not hasattr(self,"pairs"):
            self.pairs = {}
            for i,taxonA in enumerate(self.taxa):
                for j,taxonB in enumerate(self.taxa):
                    if i < j:
                        self.pairs[taxonA,taxonB] = []

                        dictA = self.get_dict(col=taxonA)
                        dictB = self.get_dict(col=taxonB)

                        for c in dictA:
                            if c in dictB:
                                valA = dictA[c]
                                valB = dictB[c]

                                for idxA in valA:
                                    for idxB in valB:
                                        dA = self[idxA,"duplicates"]
                                        dB = self[idxB,"duplicates"]
                                        if dA != 1 and dB != 1:
                                            self.pairs[taxonA,taxonB] += [(idxA,idxB)]

    def __getitem__(self,idx):
        """
        Method allows quick access to the data by passing the integer key.

        In contrast to the basic wordlist, the LexStat wordlist further allows
        to access item pairs by passing a tuple.
        """
        try:
            return self._cache[idx]
        except:
            pass

        try:
            # return full data entry as list
            out = self._data[idx]
            self._cache[idx] = out
            return out
        except:
            try:
                out = (
                        self._data[idx[0][0]][self._header[self._alias[idx[1]]]],
                        self._data[idx[0][1]][self._header[self._alias[idx[1]]]]
                        )
                return out
            except:
                try:
                    # return data entry with specified key word
                    out = self._data[idx[0]][self._header[self._alias[idx[1]]]]
                    self._cache[idx] = out
                    return out
                except:
                    pass
            
    def _get_corrdist(
            self,
            threshold = 0.7,
            modes = [("global",-2,0.5),("local",-1,0.5)],
            factor = 0.3,
            restricted_chars = '_T'
            ):
        """
        Use alignments to get a correspondences statistics.
        """
        corrdist = {}
        for i,tA in enumerate(self.taxa):
            for j,tB in enumerate(self.taxa):
                if i < j:
                    print("[i] Calculating alignments for pair {0} / {1}.".format(
                        tA,
                        tB
                        ))
                    corrdist[tA,tB] = {}
                    for mode,gop,scale in modes:
                        numbers = [self[pair,"numbers"] for pair in
                                self.pairs[tA,tB]]
                        weights = [self[pair,"weights"] for pair in
                                self.pairs[tA,tB]]
                        prostrings = [self[pair,"prostrings"] for pair in
                            self.pairs[tA,tB]]
                        corrs = calign.corrdist(
                                threshold,
                                numbers,
                                weights,
                                prostrings,
                                gop,
                                scale,
                                factor,
                                self.scorer,
                                mode,
                                restricted_chars
                                )
                                 
                        # change representation of gaps
                        for a,b in list(corrs.keys()):
                            d = corrs[a,b]
                            if a == '-':
                                a = str(i+1)+'.X.-'
                            elif b == '-':
                                b = str(j+1)+'.X.-'
                            try:
                                corrdist[tA,tB][a,b] += d / len(modes)
                            except:
                                corrdist[tA,tB][a,b] = d / len(modes)

        return corrdist

    def _get_randist(
            self,
            scaler = 5,
            threshold = 0.7,
            modes = [("global",-2,0.5),("local",-1,0.5)],
            factor = 0.3,
            restricted_chars = '_T'
            ):
        """
        Return the aligned results of randomly aligned sequences.
        """
        
        seqs = {}
        pros = {}
        weights = {}

        # determine the scale for number of random strings
        limit = self.height * scaler
        trials = self.height * scaler * 10
    
        for i,taxon in enumerate(self.taxa):
            print("[i] Analyzing taxon {0}.".format(taxon))
            tokens = self.get_list(
                    col=taxon,
                    entry="tokens",
                    flat=True
                    )
            prostrings = self.get_list(
                    col=taxon,
                    entry="prostrings",
                    flat=True
                    )
            m = MCPhon(tokens,True,prostrings)
            words = []
            j = 0
            k = 0
            while j < limit and k < trials:
                s = m.get_string()
                if s in words:
                    k += 1
                else:
                    j += 1
                    words += [s]
            
            seqs[taxon] = []
            pros[taxon] = []
            weights[taxon] = []

            for w in words:
                cls = tokens2class(w.split(' '),self.model)
                pros[taxon] += [prosodic_string(w.split(' '))]
                weights[taxon] += [prosodic_weights(pros[taxon][-1])]
                seqs[taxon] += [['{0}.{1}.{2}'.format(
                    i+1,
                    c,
                    p
                    ) for c,p in zip(cls,pros[taxon][-1])
                    ]]

        corrdist = {}
        for i,tA in enumerate(self.taxa):
            for j,tB in enumerate(self.taxa):
                if i < j:
                    print("[i] Calculating random alignments for pair {0} / {1}.".format(
                        tA,
                        tB
                        ))
                    corrdist[tA,tB] = {}
                    for mode,gop,scale in modes:
                        numbers = list(
                                zip(
                                    seqs[tA],
                                    seqs[tB]
                                        )
                                    )
                        gops = list(
                                zip(
                                    weights[tA],
                                    weights[tB]
                                    )
                                )
                        prostrings = list(
                                zip(
                                    pros[tA],
                                    pros[tB]
                                    )
                                )

                        corrs = calign.corrdist(
                                2.0,
                                numbers,
                                gops,
                                prostrings,
                                gop,
                                scale,
                                factor,
                                self.scorer,
                                mode,
                                restricted_chars
                                )
                                 
                        # change representation of gaps
                        for a,b in list(corrs.keys()):
                            d = corrs[a,b] / scaler
                            if a == '-':
                                a = str(i+1)+'.X.-'
                            elif b == '-':
                                b = str(j+1)+'.X.-'
                            try:
                                corrdist[tA,tB][a,b] += (d / len(modes)) #/ scale
                            except:
                                corrdist[tA,tB][a,b] = (d / len(modes)) #/ scale

        return corrdist








    def _align_pairs(
            self,
            mode = "global",
            gop = -2,
            scale = 0.5,
            factor = 0.3,
            restricted_chars = '_T'
            ):
        """
        Align all words for all language pairs.
        """
        # XXX note that we can improve timing here if we no longer go for
        # distances XXX
        alignments = {}
        for i,taxonA in enumerate(self.taxa):
            for j,taxonB in enumerate(self.taxa):
                if i < j:
                    # get the strings
                    numbers = [self[pair,"numbers"] for pair in
                            self.pairs[taxonA,taxonB]]
                    weights = [self[pair,"weights"] for pair in
                            self.pairs[taxonA,taxonB]]
                    prostrings = [self[pair,"prostrings"] for pair in
                        self.pairs[taxonA,taxonB]]

                    # carry out alignment
                    alms = calign.align_pairs(
                            numbers,
                            weights,
                            prostrings,
                            gop,
                            scale,
                            factor,
                            self.scorer,
                            mode,
                            restricted_chars,
                            1
                            )
                    if mode == "local":
                        alms = [(a[1],b[1],c) for a,b,c in alms]
                    alignments[taxonA,taxonB] = alms
        return alignments

