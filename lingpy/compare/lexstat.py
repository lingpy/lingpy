# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-12 11:56
# modified : 2013-04-06 16:03
"""
LexStat algorithm for automatic cognate detection.
"""

__author__="Johann-Mattis List"
__date__="2013-04-06"

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
    Basic class for automatic cognate detection.

    Parameters
    ----------
    filename : str 
        The name of the file that shall be loaded.
    model : :py:class:`~lingpy.data.model.Model` 
        The sound-class model that shall be used for the analysis. Defaults to
        the SCA sound-class model.
    merge_vowels : bool (default=True)
        Indicate whether consecutive vowels should be merged into single tokens or kept
        apart as separate tokens.
    transform : dict
        A dictionary that indicates how prosodic strings should be simplified
        (or generally transformed), using a simple key-value structure with the
        key referring to the original prosodic context and the value to the new
        value.
        Currently, prosodic strings (see
        :py:meth:`~lingpy.sequence.sound_classes.prosodic_string`) offer 11
        different prosodic contexts. Since not all these are helpful in
        preliminary analyses for cognate detection, it is useful to merge some
        of these contexts into one. The default settings distinguish only 5
        instead of 11 available contexts, namely:

        * ``C`` for all consonants in prosodically ascending position,
        * ``c`` for all consonants in prosodically descending position, 
        * ``V`` for all vowels,
        * ``T`` for all tones, and 
        * ``_`` for word-breaks.
    check : bool (default=False)
        If set to c{True}, the input file will first be checked for errors
        before the calculation is carried out. Errors will be written to the
        file ``errors.log``.

    Notes
    -----
    Instantiating this class does not require a lot of parameters. However,
    the user may modify its behaviour by providing additional attributes in the
    input file.

    """
    
    def __init__(
            self,
            filename,
            **keywords
            ):

        defaults = {
                "model" : sca,
                "merge_vowels" : True,
                'transform' : {                    
                    'A':'C', 
                    'B':'C',
                    'C':'C',
                    'L':'c',
                    'M':'c',
                    'N':'c',
                    'X':'V', #
                    'Y':'V', #
                    'Z':'V', #
                    'T':'T', #
                    '_':'_'
                    },
                "check" : False
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
            self.add_entries(
                    "tokens",
                    "ipa",
                    lambda x:ipa2tokens(
                        x,
                        merge_vowels = keywords['merge_vowels']
                        )
                    )

        # add a debug procedure for tokens
        if keywords["check"]:
            errors = []
            for key in self:
                line = self[key,"tokens"]
                if "" in line:
                    errors += [(
                        key,
                        "empty token",
                        ' '.join(line)
                        )]
                else:
                    try:
                        sonars = tokens2class(line,art)
                        if not sonars or sonars == ['0']:
                            errors += [(
                                key,
                                "empty sound-class string",
                                ' '.join(line)
                                )]
                    except:
                        errors += [(
                            key,
                            "sound-class conversion failed",
                            ' '.join(line)
                            )]
            if errors:
                out = open("errors.log","w")
                out.write("ID\tTokens\tError-Type\n")
                for a,b,c in errors:
                    out.write("{0}\t<{1}>\t{2}\n".format(a,c,b))
                out.close()
                answer = input("[?] There were errors in the input data. "
                        "Do you want to exclude the errors? (Y/N) ")
                if answer in ['Y','y','yes','j','J']:
                    self.output(
                            'csv',
                            filename=self.filename+'_cleaned',
                            subset=True,
                            rows = {"ID":"not in "+str([i[0] for i in errors])}
                            )
                    # load the data in another wordlist and copy the stuff
                    wl = Wordlist(self.filename+'_cleaned.csv')
                    
                    # change the attributes
                    self._array = wl._array
                    self._data = wl._data
                    self._dict = wl._dict
                    self._idx = wl._idx

                    # store errors in meta
                    self._meta['errors'] = [i[0] for i in errors]

                else:
                    return
            else:
                print("[i] No obvious errors found in dataset.")
        
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
            # tuples, note that this is still wip, we have to tweak around with
            # this in order to find an optimum for the calculation
            self._transform =  keywords['transform']
            #{
            #        'A':'B', 
            #        'B':'B',
            #        'C':'B',
            #        'L':'L',
            #        'M':'L',
            #        'N':'L',
            #        'X':'X', #
            #        'Y':'X', #
            #        'Z':'X', #
            #        'T':'T', #
            #        '_':'_'
            #        }
            self.add_entries(
                    "numbers",
                    "langid,classes,prostrings",
                    lambda x,y: ["{0}.{1}.{2}".format(
                        x[y[0]],
                        a,
                        self._transform[b]
                        ) for a,b in zip(x[y[1]],x[y[2]])]    
                    )

        # check for weights
        if not "weights" in self.header:
            self.add_entries(
                    "weights",
                    "prostrings",
                    lambda x:prosodic_weights(x)
                    )

        # check for duplicates
        # first, check for item 'words' in data, if this is not given, create
        # it
        if 'ipa' in self.header:
            pass
        else:
            self.add_entries('ipa','tokens',lambda x:''.join(x))

        if not "duplicates" in self.header:
            duplicates = {}
            for taxon in self.taxa:
                words = []
                for idx in self.get_list(
                        col=taxon,
                        flat=True
                        ):
                    # get the words
                    word = self[idx,'ipa']
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
            restricted_chars = '_T',
            preprocessing = False,
            gop = -2,
            cluster_method = "ugpma",
            verbose = False
            ):
        """
        Use alignments to get a correspondences statistics.
        """
        self._included = {}
        corrdist = {}

        if preprocessing:
            if 'scaid' in self.header:
                pass
            else:
                self.cluster(
                        method='sca',
                        threshold=threshold,
                        gop = gop,
                        cluster_method=cluster_method
                        )
        
        for i,tA in enumerate(self.taxa):
            for j,tB in enumerate(self.taxa):
                if i <= j:
                    if verbose: print("[i] Calculating alignments for pair {0} / {1}.".format(
                        tA,
                        tB
                        ))
                    corrdist[tA,tB] = {}
                    for mode,gop,scale in modes:
                        if preprocessing: 
                            numbers = [self[pair,"numbers"] for pair in
                                    self.pairs[tA,tB] if self[pair,"scaid"][0] == self[pair,'scaid'][1]]
                            weights = [self[pair,"weights"] for pair in
                                    self.pairs[tA,tB] if self[pair,"scaid"][0] == self[pair,'scaid'][1]]
                            prostrings = [self[pair,"prostrings"] for pair in
                                    self.pairs[tA,tB] if self[pair,"scaid"][0] == self[pair,'scaid'][1]]
                            corrs,included = calign.corrdist(
                                    10.0,
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

                        else:    
                            numbers = [self[pair,"numbers"] for pair in
                                    self.pairs[tA,tB]]
                            weights = [self[pair,"weights"] for pair in
                                    self.pairs[tA,tB]]
                            prostrings = [self[pair,"prostrings"] for pair in
                                self.pairs[tA,tB]]
                            corrs,included = calign.corrdist(
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

                        self._included[tA,tB] = included

                        # change representation of gaps
                        for a,b in list(corrs.keys()):
                            # XXX check for bias XXX
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
            method = 'markov',
            runs = 1000,
            modes = [("global",-2,0.5),("local",-1,0.5)],
            factor = 0.3,
            restricted_chars = '_T',
            rands = 1000,
            limit = 10000,
            verbose = False
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

            # get a random distribution for all pairs
            sample = random.sample(
                    [(i,j) for i in range(rands) for j in range(rands)],
                    runs
                    )

            for i,taxon in enumerate(self.taxa):
                if verbose: print("[i] Analyzing taxon {0}.".format(taxon))
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
                while j < rands:
                    s = m.get_string(new=False)
                    if s in words:
                        k += 1
                    elif k < limit:
                        j += 1
                        words += [s]
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
                        ) for c,p in zip(
                            cls,
                            [self._transform[pr] for pr in pros[taxon][-1]]
                            )
                        ]]
            
            corrdist = {}
            for i,tA in enumerate(self.taxa):
                for j,tB in enumerate(self.taxa):
                    if i <= j:
                        if verbose: print("[i] Calculating random alignments for pair {0} / {1}.".format(
                            tA,
                            tB
                            ))
                        corrdist[tA,tB] = {}
                        for mode,gop,scale in modes:
                            numbers = [(seqs[tA][x],seqs[tB][y]) for x,y in sample]
                            gops = [(weights[tA][x],weights[tB][y]) for x,y in sample]
                            prostrings = [(pros[tA][x],pros[tB][y]) for x,y in sample]

                            corrs,included = calign.corrdist(
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

                                # get the correspondence count
                                d = corrs[a,b] * self._included[tA,tB] / included # XXX check XXX * len(self.pairs[tA,tB]) / runs

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

        # use shuffle approach otherwise
        else:
            corrdist = {}
            for i,tA in enumerate(self.taxa):
                for j,tB in enumerate(self.taxa):
                    if i <= j:
                        if verbose: print("[i] Calculating random alignments for pair {0} / {1}.".format(
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

                        try:
                            sample = random.sample(
                                    [(x,y) for x in range(len(numbers)) for y in
                                        range(len(numbers))],
                                    runs
                                    )
                        # handle exception of sample is larger than population
                        except ValueError:
                            sample = [(x,y) for x in range(len(numbers)) for y
                                    in range(len(numbers))]
                        
                        # get an index that will be repeatedly changed
                        #indices = list(range(len(numbers)))

                        for mode,gop,scale in modes:
                            nnums = [(numbers[s[0]][0],numbers[s[1]][1]) for
                                    s in sample]
                            ggops = [(gops[s[0]][0],gops[s[1]][1]) for s in
                                    sample]
                            ppros = [(prostrings[s[0]][0],prostrings[s[1]][1]) for s in
                                    sample]

                            corrs,included = calign.corrdist(
                                    10.0,
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

                                # get the correspondence count
                                d = corrs[a,b] * self._included[tA,tB] / included #XXX check XXX* len(self.pairs[tA,tB]) / runs

                                # check for gaps
                                if a == '-':
                                    a = str(i+1)+'.X.-'

                                elif b == '-':
                                    b = str(j+1)+'.X.-'
                                
                                # append to overall dist
                                try:
                                    corrdist[tA,tB][a,b] += d / len(modes)
                                except:
                                    corrdist[tA,tB][a,b] = d / len(modes)
        return corrdist

    def get_scorer(
            self,
            method = 'shuffle',
            ratio = (3,2),
            vscale = 0.5,
            runs = 1000,
            threshold = 0.7,
            modes = [("global",-2,0.5),("local",-1,0.5)],
            factor = 0.3,
            restricted_chars = '_T',
            force = False,
            preprocessing = True,
            rands = 1000,
            limit = 10000,
            verbose = False,
            cluster_method = "upgma",
            gop = -2
            ):
        """
        Create a scoring function based on sound correspondences.

        Parameters
        ----------
        method : str (default='markov')
            Select between "markov", for automatically generated random
            strings, and "shuffle", for random strings taken directly from the
            data.
        ratio : tuple (default=3,2)
            Define the ratio between derived and original score for
            sound-matches.
        vscale : float (default=0.5)
            Define a scaling factor for vowels, in order to decrease their
            score in the calculations.
        runs : int (default=1000)
            Choose the number of random runs that shall be made in order to
            derive the random distribution.
        threshold : float (default=0.7)
            The threshold which used to select those words that are compared in
            order to derive the attested distribution. 
        modes : list (default = [("global",-2,0.5),("local",-1,0.5)])
            The modes which are used in order to derive the distributions from
            pairwise alignments.
        factor : float (default=0.3)
            The scaling factor for sound segments with identical prosodic
            environment.
        force : bool (default=False)
            Force recalculation of existing distribution.
        preprocessing: bool (default=False)
            Select whether SCA-analysis shall be used to derive a preliminary
            set of cognates from which the attested distribution shall be
            derived.
        rands : int (default=1000)
            If "method" is set to "markov", this parameter defines the number
            of strings to produce for the calculation of the random
            distribution.
        limit : int (default=10000)
            If "method" is set to "markov", this parameter defines the limit
            above which no more search for unique strings will be carried out.
        cluster_method : {"upgma" "single" "complete"} (default="upgma")
            Select the method to be used for the calculation of cognates in the
            preprocessing phase, if "preprocessing" is set to c{True}.
        gop : int (default=-2)
            If "preprocessing" is selected, define the gap opening penalty for
            the preprocessing calculation of cognates.
        """
        # get parameters and store them in string
        modestring = []
        for a,b,c in modes:
            modestring += ['{0}-{1}-{2:.2f}'.format(a,abs(b),c)]
        modestring = ':'.join(modestring)

        params = '{0[0]}:{0[1]}_{1:.2f}_{2}_{3:.2f}_{4}_{5:.2f}_{6}_{7}_{8}'.format(
                ratio,
                vscale,
                runs,
                threshold,
                modestring,
                factor,
                restricted_chars,
                method,
                '{0}:{1}:{2}'.format(
                    preprocessing,
                    cluster_method,
                    gop
                    )
                )

        # check for attribute
        if hasattr(self,'params') and not force:
            if self.params['scorer'] == params:
                print(
                        "[i] An identical scoring function has already been calculated, ",
                        end=''
                        )
                print("force recalculation by setting 'force' to 'True'.")
                return
            else:
                print(
                        "[i] A different scoring function has already been calculated, ",
                        end=''
                        )
                print("overwriting previous settings.") 

        # store parameters
        self.params = {'scorer':params }
        self._stamp += "# Parameters: "+params+'\n'

        # get the correspondence distribution
        corrdist = self._get_corrdist(
                threshold,
                modes,
                factor,
                restricted_chars,
                preprocessing,
                gop,
                cluster_method,
                verbose = verbose
                )
        # get the random distribution
        randist = self._get_randist(
                method,
                runs,
                modes,
                factor,
                restricted_chars,
                rands,
                limit,
                verbose = verbose
                )
        
        # store the distributions as attributes
        self._corrdist = corrdist
        self._randist = randist
        
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
                                score = np.log2((att ** 2 ) / 0.00001 )
                            elif exp and not att:
                                score = -5  #XXX gop ??? 
                            elif not exp and not att:
                                score = -90 # ???

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
                                if charA[4] in 'XYZT_' and charB[4] in 'XYZT_':
                                    matrix[idxA][idxB] = vscale * rscore
                                    matrix[idxB][idxA] = vscale * rscore
                                else:
                                    matrix[idxA][idxB] = rscore
                                    matrix[idxB][idxA] = rscore
                            except:
                                pass
        
        self.cscorer = misc.ScoreDict(self.chars,matrix)

    def align_pairs(
            self,
            idxA,
            idxB,
            method = 'lexstat',
            mode = "global",
            gop = -2,
            scale = 0.5,
            factor = 0.3,
            restricted_chars = '_T',
            distance = True,
            concept = None,
            pprint = True,
            return_distance = False,
            **keywords
            ):
        """
        Align all or some words of a given pair of languages.

        Paramters
        ---------
        idxA,idxB : {int, str}
            Use an integer to refer to the words by their unique internal ID,
            use language names to select all words for a given language.
        method : {'lexstat','sca'}
            Define the method to be used for the alignment of the words.
        mode : {'global','local','overlap','dialign'} (default='overlap')
            Select the mode for the alignment analysis.
        gop : int (default=-2)
            If 'sca' is selected as a method, define the gap opening penalty.
        scale : float (default=0.5)
            Select the scale for the gap extension penalty.
        factor : float (default=0.3)
            Select the factor for extra scores for identical prosodic segments.
        restricted_chars : str (default="T_")
            Select the restricted chars (boundary markers) in the prosodic
            strings in order to enable secondary alignment.
        distance : bool (default=True)
            If set to c{True}, return the distance instead of the similarity
            score.
        pprint : bool (default=True)
            If set to c{True}, print the results to the terminal.
        return_distance : bool (default=False)
            If set to c{True}, return the distance score, otherwise, nothing
            will be returned.
        """
        # add keywords to keywords
        keywords['method'] = method
        keywords['mode'] = mode
        keywords['scale'] = scale
        keywords['factor'] = factor
        keywords['distance'] = distance
        keywords['restricted_chars'] = restricted_chars

        if type(idxA) in [tuple,str]:
            if type(idxA) == tuple:
                idxsA = self.get_dict(col=idxA[0])[idxA[1]]
                idxsB = self.get_dict(col=idxB[0])[idxB[1]]
                for i,idxA in enumerate(idxsA):
                    for j,idxB in enumerate(idxsB):
                        self.align_pairs(idxA,idxB,**keywords)

            else:
                if not concept:
                    for c in self.concepts:
                        print("Concept: {0}".format(c))
                        keywords['concept'] = c
                        self.align_pairs(idxA,idxB,**keywords)
                        print('')
                else:
                    self.align_pairs(
                            (idxA,concept),
                            (idxB,concept),
                            concept=None,
                            **keywords
                            )
            return

        # assign the distance value
        distance = 1 if distance else 0

        # get the language ids
        lA = self[idxA,'langid']
        lB = self[idxB,'langid']

        if method == 'lexstat':
            scorer = self.cscorer
            gop = 1.0
            weightsA = [self.cscorer[str(lA)+'.X.-',n] for n in
                self[idxA,'numbers']]
            weightsB = [self.cscorer[str(lB)+'.X.-',n] for n in
                self[idxB,'numbers']]

        else:
            weightsA = self[idxA,'weights']
            weightsB = self[idxB,'weights']
            scorer = self.scorer

        almA,almB,d = calign.align_pair(
                self[idxA,'numbers'],
                self[idxB,'numbers'],
                weightsA,
                weightsB,
                self[idxA,'prostrings'],
                self[idxB,'prostrings'],
                gop,
                scale,
                factor,
                scorer,
                mode,
                restricted_chars,
                distance
                )

        # get a string of scores
        if method == 'lexstat':
            fun = lambda x,y: x if x != '-' else '{0}.X.-'.format(y)

            scoreA = [fun(a,lA) for a in almA]
            scoreB = [fun(b,lB) for b in almB]
        else:
            scoreA = almA
            scoreB = almB

        scores = ['{0:.2f}'.format(scorer[a,b]) for a,b in zip(scoreA,scoreB)]

        almA = class2tokens(self[idxA,'tokens'],almA)
        almB = class2tokens(self[idxB,'tokens'],almB)
        if pprint:
            print('\t'.join(almA))
            print('\t'.join(almB))
            print('\t'.join(scores))
            if distance:
                print('Distance: {0:.2f}'.format(d))
            else:
                print('Similarity: {0:.2f}'.format(d))
        
        if return_distance:
            return d
            
    def cluster(
            self,
            method = 'sca',
            cluster_method='upgma',
            threshold = 0.55,
            scale = 0.5,
            factor = 0.3,
            restricted_chars = '_T',
            mode = 'overlap',
            verbose = False,
            gop = -2,
            restriction = '',
            **keywords
            ):
        """
        Function for flat clustering of words into cognate sets.

        Parameters
        ----------
        method : {'sca','lexstat','edit-dist','turchin'} (default='sca')
            Select the method that shall be used for the calculation.
        cluster_method : {'upgma','single','complete'} (default='upgma')
            Select the cluster method. 'upgma' (:evobib:`Sokal1958` refers to
            average linkage clustering.
        threshold : float (default=0.6)
            Select the threshold for the cluster approach. If set to c{False},
            an automatic threshold will be calculated by calculating the
            average distance of unrelated sequences (use with care).
        scale : float (default=0.5)
            Select the scale for the gap extension penalty.
        factor : float (default=0.3)
            Select the factor for extra scores for identical prosodic segments.
        restricted_chars : str (default="T_")
            Select the restricted chars (boundary markers) in the prosodic
            strings in order to enable secondary alignment.
        mode : {'global','local','overlap','dialign'} (default='overlap')
            Select the mode for the alignment analysis.
        verbose : bool (default=False)
            Define whether verbose output should be used or not.
        gop : int (default=-2)
            If 'sca' is selected as a method, define the gap opening penalty.
        restriction : {'cv'} (default="")
            Specify the restriction for calculations using the edit-distance.
            Currently, only "cv" is supported. If *edit-dist* is selected as
            *method* and *restriction* is set to *cv*, consonant-vowel matches
            will be prohibited in the calculations and the edit distance will
            be normalized by the length of the alignment rather than the length
            of the longest sequence, as described in :evobib:`Heeringa2006`.

        """
        if not threshold:
            # use the 5 percentile of the random distribution of non-related
            # words (cross-semantic alignments) in order to determine a
            # suitable threshold for the analysis.
            if verbose: print("[i] Calculating a threshold for the calculation.")
            d = self.get_random_distances(
                    method=method,
                    mode=mode
                    )
            threshold = d[int(len(d) * 15 / 1000)]

        if hasattr(self,'params'):
            pass
        else:
            self.params = {}
        
        self.params['cluster'] = "{0}_{1}_{2:.2f}".format(
                method,
                cluster_method,
                threshold
                )
        self._stamp += '# Cluster: ' + self.params['cluster']
        
        if method not in ['lexstat','sca','turchin','edit-dist']:
            raise ValueError(
                    "[!] The method you selected is not available."
                    )

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
                entry = 'tokens'

            # define function with lamda
            function = lambda idxA,idxB: edit_dist(
                    self[idxA,entry],
                    self[idxB,entry],
                    True,
                    restriction
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

        for concept in sorted(concepts):
            if verbose: print("[i] Analyzing concept <{0}>.".format(concept))

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
            c = cluster.flat_cluster(cluster_method,threshold,matrix,revert=True)

            # extract the clusters
            clusters = [c[i]+k for i in range(len(matrix))]

            # reassign the "k" value
            k = max(clusters)
            
            # add values to cluster dictionary
            for idxA,idxB in zip(indices,clusters):
                clr[idxA] = idxB
        
        if 'override' in keywords:
            override = keywords['override']
        else:
            override = False

        if method == 'turchin':
            self.add_entries('turchinid',clr,lambda x:x,override=override)
        elif method == 'lexstat':
            self.add_entries('lexstatid',clr,lambda x:x,override=override)
        elif method == 'sca':
            self.add_entries('scaid',clr,lambda x:x,override=override)
        else:
            self.add_entries('editid',clr,lambda x:x,override=override)       
        
    def get_random_distances(
            self,
            method='lexstat',
            runs = 100,
            mode = 'overlap',
            gop = -2,
            scale = 0.5,
            factor = 0.3,
            restricted_chars = 'T_'
            ):
        """
        Method calculates randoms scores for unrelated words in a dataset.

        Parameters
        ----------
        method : {'sca','lexstat','edit-dist','turchin'} (default='sca')
            Select the method that shall be used for the calculation.
        runs : int (default=100)
            Select the number of random alignments for each language pair.
        mode : {'global','local','overlap','dialign'} (default='overlap')
            Select the mode for the alignment analysis.
        gop : int (default=-2)
            If 'sca' is selected as a method, define the gap opening penalty.
        scale : float (default=0.5)
            Select the scale for the gap extension penalty.
        factor : float (default=0.3)
            Select the factor for extra scores for identical prosodic segments.
        restricted_chars : str (default="T_")
            Select the restricted chars (boundary markers) in the prosodic
            strings in order to enable secondary alignment.

        Returns
        -------
        D : c{numpy.array}
            An array with all distances calculated for each sequence pair.
        """
        D = []
        
        if method in ['sca','lexstat']:
            function = lambda x,y: self.align_pairs(
                    x,
                    y,
                    method=method,
                    distance=True,
                    return_distance=True,
                    pprint=False,
                    mode = mode,
                    scale = scale,
                    factor = factor,
                    gop = gop
                    )
        else:
            function = lambda x,y: edit_dist(
                    self[x,'tokens'],
                    self[y,'tokens']
                    )

        for i,taxA in enumerate(self.taxa):
            for j,taxB in enumerate(self.taxa):
                if i < j:

                    # get a random selection of words from both taxa
                    pairs = self.pairs[taxA,taxB]
                    
                    try:
                        sample = random.sample(
                                [(x,y) for x in range(len(pairs)) for y in
                                    range(len(pairs))],
                                runs
                                )
                    except ValueError:
                        sample = random.sample(
                                [(x,y) for x in range(len(pairs)) for y in
                                    range(len(pairs))],
                                len(pairs)
                                )

                    sample_pairs = [(pairs[x][0],pairs[y][1]) for x,y in sample]
                    for pA,pB in sample_pairs:
                        d = function(pA,pB)

                        D += [d]

        return sorted(D)

    #def output(
    #        fileformat,
    #        **keywords
    #        ):
    #    """
    #    Write wordlist to file.

    #    Parameters
    #    ----------
    #    fileformat : {'csv', 'tre','nwk','dst', 'taxa', 'starling', 'paps.nex', 'paps.csv'}
    #        The format that is written to file. This corresponds to the file
    #        extension, thus 'csv' creates a file in csv-format, 'dst' creates
    #        a file in Phylip-distance format, etc.
    #    filename : str
    #        Specify the name of the output file (defaults to a filename that
    #        indicates the creation date).
    #    subset : bool (default=False)
    #        If set to c{True}, return only a subset of the data. Which subset
    #        is specified in the keywords 'cols' and 'rows'.
    #    cols : list
    #        If *subset* is set to c{True}, specify the columns that shall be
    #        written to the csv-file.
    #    rows : dict
    #        If *subset* is set to c{True}, use a dictionary consisting of keys
    #        that specify a column and values that give a Python-statement in
    #        raw text, such as, e.g., "== 'hand'". The content of the specified
    #        column will then be checked against statement passed in the
    #        dictionary, and if it is evaluated to c{True}, the respective row
    #        will be written to file.
    #    cognates : str
    #        Name of the column that contains the cognate IDs if 'starling' is
    #        chosen as an output format.

    #    missing : { str, int } (default=0)
    #        If 'paps.nex' or 'paps.csv' is chosen as fileformat, this character
    #        will be inserted as an indicator of missing data.

    #    tree_calc : {'neighbor', 'upgma'}
    #        If no tree has been calculated and 'tre' or 'nwk' is chosen as
    #        output format, the method that is used to calculate the tree.

    #    threshold : float (default=0.6)
    #        The threshold that is used to carry out a flat cluster analysis if
    #        'groups' or 'cluster' is chosen as output format.
    #    
    #    """
    #    
    #    if fileformat not in ['alm']:
    #        return self._output(fileformat,**keywords)
    #    
    #    if fileformat == 'alm':
    #        pass

    #        ## check for "cognates"-keywords
    #        #if "cognates" not in keywords:
    #        #    cognates = 'cogid'
    #        #for concept in self.concepts:
    #        #    l = self.get_list(
    #        #            row=concept,
    #        #            entry=cognates,
    #        #            flat=True
    #        #            )
    #        #    l = sorted(set(l))
                
