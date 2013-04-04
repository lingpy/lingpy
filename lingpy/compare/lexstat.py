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

# thirdparty
import numpy as np

# thirdparty modules
from ..thirdparty import cogent as cg

# lingpy-modules
from ..data import *
from ..sequence.sound_classes import *
from ..sequence.generate import MCPhon
from ..basic import Wordlist
from ..align.pairwise import turchin,edit_dist

try:
    from ..algorithm.cython import calign
    from ..algorithm.cython import misc
    from ..algorithm.cython import cluster
except:
    from ..algorithm.cython import _calign as calign
    from ..algorithm.cython import _misc as misc
    from ..algorithm.cython import _cluster as cluster

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

        # set the lexstat stamp
        self._stamp = "# Created using the LexStat class of LingPy-2.0\n"

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
            # change the discriminative potential of the sound-class string
            # tuples
            transform = {
                    'A':'A',
                    'B':'B',
                    'C':'C',
                    'L':'L',
                    'M':'M',
                    'N':'N',
                    'X':'X',
                    'Y':'Y',
                    'Z':'Z',
                    'T':'T',
                    '_':'_'
                    }
            self.add_entries(
                    "numbers",
                    "langid,classes,prostrings",
                    lambda x,y: ["{0}.{1}.{2}".format(x[y[0]],a,transform[b]) for a,b in zip(x[y[1]],x[y[2]])]    
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
            self.rchars = []
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
                        self.rchars.append(char[char.index('.')+1:])
            self.chars = list(set(self.chars))
            self.rchars = list(set(self.rchars))
            for i in range(self.width):
                self.chars += [str(i+1)+'.X.-']

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
        
            self.scorer = misc.ScoreDict(self.chars,matrix)

            matrix = [[0.0 for i in range(len(self.rchars))] for j in
                    range(len(self.rchars))]
            for i,charA in enumerate(self.rchars):
                for j,charB in enumerate(self.rchars):
                    if i < j:
                        
                        # add dictionary scores to the scoredict
                        score = keywords["model"](
                                charA[0],
                                charB[0]
                                )
                        matrix[i][j] = score
                        matrix[j][i] = score
                    elif i == j:
                        # add dictionary scores to the scoredict
                        score = keywords["model"](
                                charA[0],
                                charB[0]
                                )
                        matrix[i][j] = score
        
            self.rscorer = misc.ScoreDict(self.rchars,matrix)


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
                    elif i == j:
                        self.pairs[taxonA,taxonA] = []
                        dictAB = self.get_dict(col=taxonA)
                        for c in dictAB:
                            valAB = dictAB[c]
                            for idx in valAB:
                                dAB = self[idx,"duplicates"]
                                if dAB != 1:
                                    self.pairs[taxonA,taxonA] += [(idx,idx)]

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
                elif i == j:
                    corrdist[tA,tB] = {}

                    numbers = [self[pair,"numbers"] for pair in
                            self.pairs[tA,tA]]
                    for a,b in numbers:
                        l = len(a)
                        for k in range(l):
                            d = self.scorer[a[k],b[k]]
                            try:
                                corrdist[tA,tB][a[k],b[k]] += d / len(modes)
                            except:
                                corrdist[tA,tB][a[k],b[k]] = d / len(modes)

        return corrdist

    def _get_randist(
            self,
            method = 'markov',
            runs = 1000,
            modes = [("global",-2,0.5),("local",-1,0.5)],
            factor = 0.3,
            restricted_chars = '_T'
            ):
        """
        Return the aligned results of randomly aligned sequences.
        """
        # determine the mode
        if method in ['markov','markov-chain','mc']:
            method = 'markov'
        else:
            method = 'shuffle'
        
        if method == 'markov':
            seqs = {}
            pros = {}
            weights = {}
            
            # determine the scale for number of random strings
            #limit = self.height * scaler
            trials = runs * 5

            # get a random distribution for all pairs

            sample = random.sample(
                    [(i,j) for i in range(1000) for j in range(1000)],
                    runs
                    )

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
                while j < 1000: # and k < 10000:
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
                    seqs[taxon] += [['{0}.{1}'.format(
                        c,
                        p
                        ) for c,p in zip(cls,pros[taxon][-1])
                        ]]

            corrdist = {}
            for i,tA in enumerate(self.taxa):
                for j,tB in enumerate(self.taxa):
                    if i <= j:
                        print("[i] Calculating random alignments for pair {0} / {1}.".format(
                            tA,
                            tB
                            ))
                        corrdist[tA,tB] = {}
                        for mode,gop,scale in modes:
                            #numbers = list(
                            #        zip(
                            #            seqs[tA],
                            #            seqs[tB]
                            #                )
                            #            )
                            #gops = list(
                            #        zip(
                            #            weights[tA],
                            #            weights[tB]
                            #            )
                            #        )
                            #prostrings = list(
                            #        zip(
                            #            pros[tA],
                            #            pros[tB]
                            #            )
                            #        )
                            numbers = [(seqs[tA][x],seqs[tB][y]) for x,y in sample]
                            gops = [(weights[tA][x],weights[tB][y]) for x,y in sample]
                            prostrings = [(pros[tA][x],pros[tB][y]) for x,y in sample]

                            corrs = calign.corrdist(
                                    10.0,
                                    numbers,
                                    gops,
                                    prostrings,
                                    gop,
                                    scale,
                                    factor,
                                    self.rscorer,
                                    mode,
                                    restricted_chars
                                    )
                                     
                            # change representation of gaps
                            for a,b in list(corrs.keys()):

                                # get the corresondence count
                                d = corrs[a,b] * len(self.pairs[tA,tB]) / runs

                                # check for gaps
                                if a == '-':
                                    a = 'X.-'
                                elif b == '-':
                                    b = 'X.-'
                                
                                a = str(i+1)+'.'+a
                                b = str(j+1)+'.'+b

                                # append to overall dist
                                try:
                                    corrdist[tA,tB][a,b] += d / len(modes)
                                except:
                                    corrdist[tA,tB][a,b] = d / len(modes)
        else:
            corrdist = {}
            for i,tA in enumerate(self.taxa):
                for j,tB in enumerate(self.taxa):
                    if i <= j:
                        print("[i] Calculating random alignments for pair {0} / {1}.".format(
                            tA,
                            tB
                            ))
                        corrdist[tA,tB] = {}

                        # get the number pairs etc.
                        numbers = [self[pair,'numbers'] for pair in
                                self.pairs[tA,tB]]
                        gops = [self[pair,'weights'] for pair in
                                self.pairs[tA,tB]]
                        prostrings = [self[pair,'prostrings'] for pair in
                                self.pairs[tA,tB]]
                        
                        # get an index that will be repeatedly changed
                        indices = list(range(len(numbers)))

                        for mode,gop,scale in modes:
                            for k in range(runs):
                                random.shuffle(indices)
                                nnums = [(numbers[indices[x]][0],numbers[x][1]) for x in
                                    range(len(indices))]
                                ggops = [(gops[indices[x]][0],gops[x][1]) for x in
                                    range(len(indices))]
                                ppros = [(prostrings[indices[x]][0],prostrings[x][1]) for x in
                                    range(len(indices))]
                                corrs = calign.corrdist(
                                        2.0,
                                        nnums,
                                        ggops,
                                        ppros,
                                        gop,
                                        scale,
                                        factor,
                                        self.scorer,
                                        mode,
                                        restricted_chars
                                        )
                                         
                                # change representation of gaps
                                for a,b in list(corrs.keys()):

                                    # get the corresondence count
                                    d = corrs[a,b] / runs #* len(self.pairs[tA,tB]) / len(numbers)

                                    # check for gaps
                                    if a == '-':
                                        a = str(i)+'.X.-'

                                    elif b == '-':
                                        b = str(j)+'X.-'
                                    
                                    # append to overall dist
                                    try:
                                        corrdist[tA,tB][a,b] += d / len(modes)
                                    except:
                                        corrdist[tA,tB][a,b] = d / len(modes)
        return corrdist

    def _get_scorer(
            self,
            method = 'markov',
            ratio = (3,2),
            vscale = 0.5,
            runs = 1000,
            threshold = 0.7,
            modes = [("global",-2,0.5),("local",-1,0.5)],
            factor = 0.3,
            restricted_chars = '_T',
            force = False
            ):
        """
        Create a scoring function based on sound correspondences.
        """
        # get parameters and store them in string
        modestring = []
        for a,b,c in modes:
            modestring += ['{0}-{1}-{2:.2f}'.format(a,abs(b),c)]
        modestring = ':'.join(modestring)

        params = '{0[0]}:{0[1]}_{1:.2f}_{2}_{3:.2f}_{4}_{5:.2f}_{6}_{7}'.format(
                ratio,
                vscale,
                runs,
                threshold,
                modestring,
                factor,
                restricted_chars,
                method
                )

        # check for attribute
        if hasattr(self,'params') and not force:
            if self.params == params:
                print("[i] An identical scoring function has already been calculated, force recalculation by setting 'force' to 'True'.")
                return
            else:
                print("[i] A different scoring function has already been calculated, overwriting previous settings.") 

        # store parameters
        self.params = params
        self._stamp += "# Parameters: "+params+'\n'

        # get the correspondence distribution
        corrdist = self._get_corrdist(
                threshold,
                modes,
                factor,
                restricted_chars
                )
        # get the random distribution
        randist = self._get_randist(
                method,
                runs,
                modes,
                factor,
                restricted_chars
                )
        
        # get the average gop
        gop = sum([m[1] for m in modes]) / len(modes)

        # create the new scoring matrix
        matrix = [[c for c in line] for line in self.scorer.matrix]
        char_dict = self.scorer.chars2int

        # start the calculation
        for i,tA in enumerate(self.taxa):
            for j,tB in enumerate(self.taxa):
                if i <= j:
                    for charA in list(self.freqs[tA]) + [str(i+1)+'.X.-']:
                        for charB in list(self.freqs[tB]) + [str(j+1)+'.X.-']:
                            try:
                                exp = randist[tA,tB][charA,charB]
                            except:
                                exp = False
                            try:
                                att = corrdist[tA,tB][charA,charB]
                            except:
                                att = False

                            # in the following we follow the former lexstat
                            # protocol
                            if att <= 1 and i != j:
                                att = False

                            if att and exp:
                                score = np.log2((att ** 2 ) / ( exp ** 2 ) )
                            elif att and not exp:
                                score = np.log2((att ** 2 ) / 0.01 )
                            elif exp and not att:
                                score = gop #XXX
                            elif not exp and not att:
                                score = -250.0

                            # combine the scores
                            if '-' not in charA+charB:
                                sim = self.scorer[charA,charB]
                            else:
                                sim = gop

                            # get the real score
                            rscore = ( ratio[0] * score + ratio[1] * sim ) / sum (ratio)
                            
                            try:
                                idxA = char_dict[charA]
                                idxB = char_dict[charB]

                                # use the vowel scale
                                if charA[2] in 'XYZT_' and charB[2] in 'XYZT_':
                                    matrix[idxA][idxB] = vscale * rscore
                                    matrix[idxB][idxA] = vscale * rscore
                                else:
                                    matrix[idxA][idxB] = rscore
                                    matrix[idxB][idxA] = rscore
                            except:
                                pass
        
        self.cscorer = misc.ScoreDict(self.chars,matrix)

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
                    prostrings = [self[pair,"prostrings"] for pair in
                        self.pairs[taxonA,taxonB]]

                    # note that the weights have to be calculated differently,
                    # once alignments scores are already calculated
                    weights = []
                    for nA,nB in numbers:
                        weightA = [self.cscorer[str(j+1)+'.X.-',n] for n in nA]
                        weightB = [self.cscorer[str(i+1)+'.X.-',n] for n in nB]
                        weights += [(weightA,weightB)]
  
                    # carry out alignment
                    alms = calign.align_pairs(
                            numbers,
                            weights,
                            prostrings,
                            1.0,
                            scale,
                            factor,
                            self.cscorer,
                            mode,
                            restricted_chars,
                            2
                            )
                    if mode == "local":
                        alms = [(a[1],b[1],c,d) for a,b,c,d in alms]
                    alignments[taxonA,taxonB] = alms
        return alignments
    
    def cluster(
            self,
            method = 'sca',
            threshold = 0.55,
            scale = 0.5,
            factor = 0.3,
            restricted_chars = '_T',
            mode = 'overlap',
            verbose = False,
            gop = -2,
            **keywords
            ):
        """
        Internal function for clustering using the LexStat approach.
        """
        # check for method
        if method == 'lexstat':
            
            # check for scorer
            if not hasattr(self,'cscorer'):
                print("[i] No correspondence-scorer has been specied.")
                return
            
            # define the function with help of lambda
            function = lambda idxA,idxy: calign.align_pair(
                    self[idxA,'numbers'],
                    self[idxB,'numbers'],
                    [self.cscorer[self[idxB,'langid'] + ".X.-",n] for n in
                        self[idxA,'numbers']],
                    [self.cscorer[self[idxA,'langid'] + ".X.-",n] for n in
                        self[idxB,'numbers']],

                    self[idxA,'prostrings'],
                    self[idxB,'prostrings'],
                    1,
                    scale,
                    factor,
                    self.cscorer,
                    mode,
                    restricted_chars,
                    1
                    )[2]
        elif method == 'sca':
            # define the function with help of lambda
            function = lambda idxA,idxB: calign.align_pair(
                    self[idxA,'numbers'],
                    self[idxB,'numbers'],
                    self[idxA,'weights'],
                    self[idxB,'weights'],
                    self[idxA,'prostrings'],
                    self[idxB,'prostrings'],
                    gop,
                    scale,
                    factor,
                    self.scorer,
                    mode,
                    restricted_chars,
                    1
                    )[2]  

        elif method == 'edit-dist':
            try:
                entry = keywords['entry']
            except:
                entry = 'word'

            # define function with lamda
            function = lambda idxA,idxB: edit_dist(
                    self[idxA,entry],
                    self[idxB,entry],
                    True
                    )

        elif method == 'turchin':
            function = lambda idxA,idxB: turchin(
                    self[idxA,'tokens'],
                    self[idxB,'tokens']
                    )
        
        # for convenience and later addons
        concepts = self.concepts

        # make a dictionary that stores the clusters for later update
        clr = {}
        k = 0

        for concept in concepts:
            if verbose: print("[i] Analyzing concept {0}.".format(concept))

            indices = self.get_list(
                    row=concept,
                    flat=True
                    )

            matrix = []
            
            for i,idxA in enumerate(indices):
                for j,idxB in enumerate(indices):
                    if i < j:
                        d = function(idxA,idxB)
                        ## append distance score to matrix
                        matrix += [d]
            matrix = misc.squareform(matrix)
            
            # calculate the clusters using flat-upgma
            c = cluster.flat_upgma(threshold,matrix,revert=True)

            # extract the clusters
            clusters = [c[i]+k for i in range(len(matrix))]

            # reassign the "k" value
            k = max(clusters)
            
            # add values to cluster dictionary
            for idxA,idxB in zip(indices,clusters):
                clr[idxA] = idxB
        
        if method == 'turchin':
            self.add_entries('turchinid',clr,lambda x:x)
        elif method == 'lexstat':
            self.add_entries('lexstatid',clr,lambda x:x)
        elif method == 'sca':
            self.add_entries('scaid',clr,lambda x:x)
        else:
            self.add_entries('editid',clr,lambda x:x)       
        
        # return the dictionary
        #return clr

    #def _get_distance(
    #        self,
    #        idxA,
    #        idxB,
    #        method,
    #        scorer
    #        ):
    #    """
    #    Internal function returns a distance from the data, depending on the
    #    method.
    #    """
    #    
    #    numA = self[idxA,'numbers']
    #    numB = self[idxB,'numbers']
    #    proA = self[idxA,'prostrings']
    #    proB = self[idxB,'prostrings']
    #    
    #    # get the weights
    #    wA = self[idxA,'weights']
    #    wB = self[idxB,'weights']
    #    almA,almB,d = calign.align_pair(
    #            numA,
    #            numB,
    #            wA,
    #            wB,
    #            proA,
    #            proB,
    #            gop,
    #            scale,
    #            factor,
    #            scorer,
    #            mode,
    #            restricted_chars,
    #            1
    #            )
    #    return d

    #def _sca(
    #        threshold = 0.55,
    #        scale = 0.5,
    #        factor = 0.3,
    #        restricted_chars = '_T',
    #        mode = 'overlap',
    #        gop = -1.0,
    #        verbose = True,
    #        **keywords,
    #        ):
    #    """
    #    Carry out an analysis using the SCA distance.
    #    """

    #    # for convenience and later addons
    #    concepts = self.concepts

    #    # make a dictionary that stores the clusters for later update
    #    clr = {}
    #    k = 0

    #    for concept in concepts:
    #        print("[i] Analyzing concept {0}.".format(concept))

    #        indices = self.get_list(
    #                row=concept,
    #                flat=True
    #                )

    #        matrix = []
    #        
    #        for i,idxA in enumerate(indices):
    #            for j,idxB in enumerate(indices):
    #                if i < j:
    #                    numA = self[idxA,'numbers']
    #                    numB = self[idxB,'numbers']
    #                    proA = self[idxA,'prostrings']
    #                    proB = self[idxB,'prostrings']
    #                    
    #                    # get language ids
    #                    lA = self[idxA,'langid']
    #                    lB = self[idxB,'langid']
    #                    
    #                    # get the weights
    #                    wA = self[idxA,'weights']
    #                    wB = self[idxB,'weights']
    #                    almA,almB,d = calign.align_pair(
    #                            numA,
    #                            numB,
    #                            wA,
    #                            wB,
    #                            proA,
    #                            proB,
    #                            gop,
    #                            scale,
    #                            factor,
    #                            scorer,
    #                            mode,
    #                            restricted_chars,
    #                            1
    #                            )

    #                    # append distance score to matrix
    #                    matrix += [d]
    #        matrix = misc.squareform(matrix)
    #        
    #        # calculate the clusters using flat-upgma
    #        c = cluster.flat_upgma(threshold,matrix,revert=True)

    #        # extract the clusters
    #        clusters = [c[i]+k for i in range(len(matrix))]

    #        # reassign the "k" value
    #        k = max(clusters)
    #        
    #        # add values to cluster dictionary
    #        for idxA,idxB in zip(indices,clusters):
    #            clr[idxA] = idxB
    #    return clr

   
    #def _turchin(self):
    #    """
    #    Calculate distances on the basis of Turchin's et al. approach.
    #    """
    #    
    #    pass

    #def _edit_dist(
    #        threshold = 0.55,
    #        verbose = True,
    #        entry = False,
    #        **keywords,
    #        ):
    #    """
    #    Carry out an analysis using the edit distance.
    #    """

    #    # check for classes
    #    if not entry:
    #        entry = 'words'

    #    # for convenience and later addons
    #    concepts = self.concepts

    #    # make a dictionary that stores the clusters for later update
    #    clr = {}
    #    k = 0

    #    for concept in concepts:
    #        print("[i] Analyzing concept {0}.".format(concept))

    #        indices = self.get_list(
    #                row=concept,
    #                flat=True
    #                )

    #        matrix = []
    #        
    #        for i,idxA in enumerate(indices):
    #            for j,idxB in enumerate(indices):
    #                if i < j:
    #                    wordA = self[idxA,entry]
    #                    wordB = self[idxB,entry]
    #                    
    #                    # get language ids
    #                    lA = self[idxA,'langid']
    #                    lB = self[idxB,'langid']
    #                    
    #                    # get the distance
    #                    d = edit_dist(wordA,wordB,normalized=True)

    #                    # append distance score to matrix
    #                    matrix += [d]
    #        matrix = misc.squareform(matrix)
    #        
    #        # calculate the clusters using flat-upgma
    #        c = cluster.flat_upgma(threshold,matrix,revert=True)

    #        # extract the clusters
    #        clusters = [c[i]+k for i in range(len(matrix))]

    #        # reassign the "k" value
    #        k = max(clusters)
    #        
    #        # add values to cluster dictionary
    #        for idxA,idxB in zip(indices,clusters):
    #            clr[idxA] = idxB

    #    return clr
