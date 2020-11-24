# *-* coding: utf-8 *-*
"""
Module provides a class for partial cognate detection, expanding the LexStat class.

"""
from __future__ import print_function, division, unicode_literals
from collections import defaultdict
from itertools import combinations, product
import random

import numpy as np
import networkx as nx

import lingpy
from lingpy.settings import rcParams
from lingpy.algorithm import clustering, extra, misc
from lingpy.algorithm import calign
from lingpy.compare.lexstat import LexStat
from lingpy import util, log 

try:
    from lingpy.algorithm.cython import calign
except ImportError:
    from lingpy.algorithm.cython import _calign as calign

# taking functions from lexstat source code here
def _charstring(id_, char='X', cls='-'):
    return '{0}.{1}.{2}'.format(id_, char, cls)

def _get_slices(tokens, **keywords):
    """
    Function returns a slice tuple that indicates the boundaries for words with\
            morpheme annotations.
    """
    kw = dict(
            sep=lingpy.settings.rcParams['morpheme_separator'],
            word_sep=lingpy.settings.rcParams['word_separator'],
            word_seps=lingpy.settings.rcParams['word_separators'],
            seps=lingpy.settings.rcParams['morpheme_separators'],
            tones='T',
            split_on_tones=False,
            )
    kw.update(keywords)
    morphemes = lingpy.basictypes.lists(tokens).n if not kw['split_on_tones'] \
            else lingpy.sequence.sound_classes.tokens2morphemes(tokens, **kw)
    out = []
    current = 0
    for morpheme in morphemes:
        out += [(current, current+len(morpheme))]
        current = current+len(morpheme)+(1 if not kw['split_on_tones'] else 0)
    return out

class Partial(LexStat):
    """
    Extended class for automatic detection of partial cognates.

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
        value.  Currently, prosodic strings (see
        :py:meth:`~lingpy.sequence.sound_classes.prosodic_string`) offer 11
        different prosodic contexts. Since not all these are helpful in
        preliminary analyses for cognate detection, it is useful to merge
        some of these contexts into one. The default settings distinguish only
        5 instead of 11 available contexts, namely:

        * ``C`` for all consonants in prosodically ascending position,
        * ``c`` for all consonants in prosodically descending position,
        * ``V`` for all vowels,
        * ``T`` for all tones, and
        * ``_`` for word-breaks.

        Make sure to check also the "vowel" keyword when initialising a LexStat
        object, since the symbols you use for vowels and tones should be
        identical with the ones you define in your transform dictionary.
    vowels : str (default="VT_")
        For scoring function creation using the
        :py:class:`~lingpy.compare.lexstat.LexStat.get_scorer` function, you
        have the possibility to use reduced scores for the matching of tones
        and vowels by modifying the "vscale" parameter, which is set to 0.5 as
        a default.  In order to make sure that vowels and tones are properly
        detected, make sure your prosodic string representation of vowels
        matches the one in this keyword. Thus, if you change the prosodic
        strings using the "transform" keyword, you also need to change the
        vowel string, to make sure that "vscale" works as wanted in the
        :py:class:`~lingpy.compare.lexstat.LexStat.get_scorer` function.
    check : bool (default=False)
        If set to **True**, the input file will first be checked for errors
        before the calculation is carried out. Errors will be written to the
        file **errors**, defaulting to ``errors.log``. See also ``apply_checks``
    apply_checks : bool (default=False)
        If set to **True**, any errors identified by `check` will be handled
        silently.
    no_bscorer: bool (default=False)
        If set to **True**, this will suppress the creation of a
        language-specific scoring function (which may become quite large and is
        additional ballast if the method "lexstat" is not used after all. If
        you use the "lexstat" method, however, this needs to be set to
        **False**.
    errors : str
        The name of the error log.

    Attributes
    ----------
    pairs : dict
        A dictionary with tuples of language names as key and indices as value, \
        pointing to unique combinations of words with the same meaning in all \
        language pairs.
    model : :py:class:`~lingpy.data.model.Model`
        The sound class model instance which serves to convert the phonetic
        data into sound classes.
    chars : list
        A list of all unique language-specific character types in the
        instantiated LexStat object. The characters in this list consist of

        * the language identifier (numeric, referenced as "langid" as a
          default, but customizable via the keyword "langid")
        * the sound class symbol for the respective IPA transcription value
        * the prosodic class value

        All values are represented in the above order as one string, separated
        by a dot. Gaps are also included in this collection. They are
        traditionally represented as "X" for the sound class and "-" for the
        prosodic string.
    rchars : list
        A list containing all unique character types across languages. In
        contrast to the chars-attribute, the "rchars" (raw chars) do not
        contain the language identifier, thus they only consist of two values,
        separated by a dot, namely, the sound class symbol, and the prosodic
        class value.
    scorer : dict
        A collection of :py:class:`~lingpy.algorithm.cython.misc.ScoreDict`
        objects, which are used to score the strings. LexStat distinguishes two
        different scoring functions:

        * rscorer: A "raw" scorer that is not language-specific and consists
          only of sound class values and prosodic string values. This scorer is
          traditionally used to carry out the first alignment in order to
          calculate the language-specific scorer. It is directly accessible as an
          attribute of the LexStat class
          (:py:class:`~lingpy.compare.lexstat.lexstat.rscorer`). The characters
          which constitute the values in this scorer are accessible via the
          "rchars" attribue of each lexstat class.
        * bscorer: The language-specific scorer. This scorer is made of unique
          language-specific characters. These are accessible via the "chars"
          attribute of each LexStat class. As the "rscorer", the "bscorer" can
          also be accessed directly as an attribute of the LexStat class 
          (:py:class:`~lingpy.compare.lexstat.lexstat.bscorer`).

    Notes
    -----
    This method automatically infers partial cognate sets from data which was
    previously morphologically segmented. 

    """

    def __init__(self, infile, **keywords):

        kw = {
                "morphemes": "morphemes",
                "partial_cognates": "partial_cognate_sets",
                "split_on_tones": False
                }
        kw.update(keywords)
        lingpy.compare.lexstat.LexStat.__init__(self, infile, **kw)
        self._morphemes = kw['morphemes']
        self._partials = kw['partial_cognates']

        # get the slices
        self._slices = {}
        for idx, tokens in self.iter_rows(self._segments):
            self._slices[idx] = _get_slices(tokens, **keywords)


    def _get_partial_corrdist(self, **keywords):
        """
        Use alignments to get a correspondences statistics.
        """
        kw = dict(
            cluster_method='upgma',
            factor=rcParams['align_factor'],
            gop=rcParams['align_gop'],
            modes=rcParams['lexstat_modes'],
            preprocessing=False,
            preprocessing_method=rcParams['lexstat_preprocessing_method'],
            preprocessing_threshold=rcParams[
                'lexstat_preprocessing_threshold'],
            split_on_tones=False,
            ref='scaid',
            restricted_chars=rcParams['restricted_chars'],
            threshold=rcParams['lexstat_scoring_threshold'],
            subset=False)
        kw.update(keywords)

        self._included = {}
        corrdist = {}

        if kw['preprocessing']:
            if kw['ref'] not in self.header:
                self.cluster(
                    method=kw['preprocessing_method'],
                    threshold=kw['preprocessing_threshold'],
                    gop=kw['gop'],
                    cluster_method=kw['cluster_method'],
                    ref=kw['ref'])

        with util.pb(
                desc='CORRESPONDENCE CALCULATION',
                total=self.width ** 2 / 2) as pb:
            for (i, tA), (j, tB) in util.multicombinations2(
                    enumerate(self.cols)):
                pb.update(1)
                log.info("Calculating alignments for pair {0} / {1}.".format(
                    tA, tB))

                corrdist[tA, tB] = defaultdict(float)
                for mode, gop, scale in kw['modes']:
                    pairs = self.pairs[tA, tB]
                    if kw['subset']:
                        pairs = [
                                pair for pair in pairs if pair in
                                self.subsets[tA, tB]]

                    # threshold and preprocessing, make sure threshold is
                    # different from pre-processing threshold when
                    # preprocessing is set to false
                    if kw['preprocessing']:
                        pairs = [pair for pair in pairs
                                 if self[pair, kw['ref']][0] == self[
                                     pair, kw['ref']][1]]
                        threshold = 10.0
                    else:
                        threshold = kw['threshold']

                    # create morpheme-segmented pairs
                    new_nums, new_weights, new_pros = [], [], []
                    for idxA, idxB in pairs:
                        for iA, iB in self._slices[idxA]:
                            for jA, jB in self._slices[idxB]:
                                new_nums += [(
                                    self[idxA, self._numbers][iA:iB],
                                    self[idxB, self._numbers][jA:jB]
                                    )]
                                new_weights += [(
                                    self[idxA, self._weights][iA:iB],
                                    self[idxB, self._weights][jA:jB]
                                    )]
                                new_pros += [(
                                    self[idxA, self._prostrings][iA:iB],
                                    self[idxB, self._prostrings][jA:jB]
                                    )]

                    corrs, self._included[tA, tB] = calign.corrdist(
                        threshold,
                        new_nums,
                        new_weights,
                        new_pros,
                        gop,
                        scale,
                        kw['factor'],
                        self.bscorer,
                        mode,
                        kw['restricted_chars'])

                    # change representation of gaps
                    for (a, b), d in corrs.items():
                        # XXX check for bias XXX
                        if a == '-':
                            a = util.charstring(i + 1)
                        elif b == '-':
                            b = util.charstring(j + 1)
                        corrdist[tA, tB][a, b] += d / float(len(kw['modes']))

        return corrdist

    def _get_partial_randist(self, **keywords):
        """
        Return the aligned results of randomly aligned sequences.
        """
        kw = dict(
            modes=rcParams['lexstat_modes'],
            factor=rcParams['align_factor'],
            restricted_chars=rcParams['restricted_chars'],
            runs=rcParams['lexstat_runs'],
            rands=rcParams['lexstat_rands'],
            limit=rcParams['lexstat_limit'],
            method=rcParams['lexstat_scoring_method'])
        kw.update(keywords)

        # determine the mode
        method = 'markov' if kw['method'] in ['markov', 'markov-chain', 'mc'] \
            else 'shuffle'

        corrdist = {}
        tasks = (self.width ** 2) / 2
        with util.pb(
                desc='RANDOM CORRESPONDENCE CALCULATION',
                total=tasks) as progress:
            for (i, tA), (j, tB) in util.multicombinations2(
                    enumerate(self.cols)):
                progress.update(1)
                log.info(
                    "Calculating random alignments"
                    "for pair {0}/{1}.".format(tA, tB)
                )
                corrdist[tA, tB] = defaultdict(float)
                
                # create morpheme-segmented pairs
                pairs = self.pairs[tA, tB]
                new_nums, new_weights, new_pros = [], [], []
                for idxA, idxB in pairs:
                    for iA, iB in self._slices[idxA]:
                        for jA, jB in self._slices[idxB]:
                            new_nums += [(
                                self[idxA, self._numbers][iA:iB],
                                self[idxB, self._numbers][jA:jB]
                                )]
                            new_weights += [(
                                self[idxA, self._weights][iA:iB],
                                self[idxB, self._weights][jA:jB]
                                )]
                            new_pros += [(
                                self[idxA, self._prostrings][iA:iB],
                                self[idxB, self._prostrings][jA:jB]
                                )]
                # get the number pairs etc.
                sample = [
                        (x, y)
                        for x in range(len(new_nums)) for y in
                        range(len(new_nums))]
                if len(sample) > kw['runs']:
                    sample = random.sample(sample, kw['runs'])

                for mode, gop, scale in kw['modes']:
                    corrs, included = calign.corrdist(
                        10.0,
                        [(
                            new_nums[s[0]][0],
                            new_nums[s[1]][1]) for s in sample],
                        [(new_weights[s[0]][0], new_weights[s[1]][1]) for s in sample],
                        [(
                            new_pros[s[0]][0],
                            new_pros[s[1]][1]) for s in sample],
                        gop,
                        scale,
                        kw['factor'],
                        self.bscorer,
                        mode,
                        kw['restricted_chars'])

                    # change representation of gaps
                    for a, b in list(corrs.keys()):
                        # get the correspondence count
                        d = corrs[a, b] * self._included[tA, tB] / included
                        # XXX check XXX* len(self.pairs[tA,tB]) / runs

                        # check for gaps
                        if a == '-':
                            a = util.charstring(i + 1)
                        elif b == '-':
                            b = util.charstring(j + 1)

                        corrdist[tA, tB][a, b] += d / len(kw['modes'])
        return corrdist

    def get_partial_scorer(self, **keywords):
        """
        Create a scoring function based on sound correspondences.

        Parameters
        ----------
        method : str (default='shuffle')
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
            The threshold which used to select those words that are compared
            in order to derive the attested distribution.
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
        unattested : {int, float} (default=-5)
            If a pair of sounds is not attested in the data, but expected by
            the alignment algorithm that computes the expected distribution,
            the score would be -infinity. Yet in order to allow to smooth this
            behaviour and to reduce the strictness, we set a default negative
            value which does not necessarily need to be too high, since it may
            well be that we miss a potentially good pairing in the first runs
            of alignment analyses. Use this keyword to adjust this parameter.
        unexpected : {int, float} (default=0.000001)
            If a pair is encountered in a given alignment but not expected
            according to the randomized alignments, the score would be not
            calculable, since we had to divide by zero. For this reason, we set
            a very small constant, by which the score is divided in this case.
            Not that this constant is only relevant in those cases where the
            shuffling procedure was not carried out long enough.

        """
        kw = dict(
            method=rcParams['lexstat_scoring_method'],
            ratio=rcParams['lexstat_ratio'],
            vscale=rcParams['lexstat_vscale'],
            runs=rcParams['lexstat_runs'],
            threshold=rcParams['lexstat_scoring_threshold'],
            modes=rcParams['lexstat_modes'],
            factor=rcParams['align_factor'],
            restricted_chars=rcParams['restricted_chars'],
            force=False,
            preprocessing=False,
            rands=rcParams['lexstat_rands'],
            limit=rcParams['lexstat_limit'],
            cluster_method=rcParams['lexstat_cluster_method'],
            gop=rcParams['align_gop'],
            preprocessing_threshold=rcParams[
                'lexstat_preprocessing_threshold'],
            preprocessing_method=rcParams['lexstat_preprocessing_method'],
            subset=False,
            defaults=False,
            unattested=-5,
            unexpected=0.00001,
            smooth=1
        )
        kw.update(keywords)
        if kw['defaults']:
            return kw

        # get parameters and store them in string
        params = dict(
            ratio=kw['ratio'],
            vscale=kw['vscale'],
            runs=kw['runs'],
            scoring_threshold=kw['threshold'],
            preprocessing_threshold=kw['preprocessing_threshold'],
            modestring=':'.join(
                '{0}-{1}-{2:.2f}'.format(a, abs(b), c) for a, b, c in
                kw['modes']),
            factor=kw['factor'],
            restricted_chars=kw['restricted_chars'],
            method=kw['method'],
            preprocessing='{0}:{1}:{2}'.format(
                kw['preprocessing'], kw['cluster_method'], kw['gop']),
            unattested=kw['unattested'],
            unexpected=kw['unexpected'],
            smooth=kw['smooth']
            )

        parstring = '_'.join(
            [
                '{ratio[0]}:{ratio[1]}',
                '{vscale:.2f}',
                '{runs}',
                '{scoring_threshold:.2f}',
                '{modestring}',
                '{factor:.2f}',
                '{restricted_chars}',
                '{method}',
                '{preprocessing}',
                '{preprocessing_threshold}',
                '{unexpected:.2f}',
                '{unattested:.2f}'
            ]).format(**params)

        # check for existing attributes
        if hasattr(self, 'cscorer') and not kw['force']:
            log.warning(
                "An identical scoring function has already been calculated, "
                "force recalculation by setting 'force' to 'True'.")
            return

        # check for attribute
        if hasattr(self, 'params') and not kw['force']:
            if 'cscorer' in self.params:
                if self.params['cscorer'] == params:
                    log.warning(
                        "An identical scoring function has already been "
                        "calculated, force recalculation by setting 'force'"
                        " to 'True'.")
                    return
            else:
                log.warning(
                    "A different scoring function has already been calculated, "
                    "overwriting previous settings.")

        # store parameters
        self.params = {'cscorer': params}
        self._meta['params'] = self.params
        self._stamp += "# Parameters: " + parstring + '\n'

        # get the correspondence distribution
        self._corrdist = self._get_partial_corrdist(**kw)
        # get the random distribution
        self._randist = self._get_partial_randist(**kw)

        # get the average gop
        gop = sum([m[1] for m in kw['modes']]) / len(kw['modes'])

        # create the new scoring matrix
        matrix = [[c for c in line] for line in self.bscorer.matrix]
        char_dict = self.bscorer.chars2int

        for (i, tA), (j, tB) in util.multicombinations2(enumerate(self.cols)):
            for charA, charB in product(
                list(self.freqs[tA]) + [util.charstring(i + 1)],
                list(self.freqs[tB]) + [util.charstring(j + 1)]
            ):
                exp = self._randist.get(
                        (tA, tB), {}).get((charA, charB), False)
                att = self._corrdist.get(
                        (tA, tB), {}).get((charA, charB), False)
                # in the following we follow the former lexstat protocol
                if att <= kw['smooth'] and i != j:
                    att = False

                if att and exp:
                    score = np.log2((att ** 2) / (exp ** 2))
                elif att and not exp:
                    score = np.log2((att ** 2) / kw['unexpected'])
                elif exp and not att:
                    score = kw['unattested']  # XXX gop ???
                else:  # elif not exp and not att:
                    score = -90  # ???

                # combine the scores
                if rcParams['gap_symbol'] not in charA + charB:
                    sim = self.bscorer[charA, charB]
                else:
                    sim = gop

                # get the real score
                rscore = (kw['ratio'][0] * score + kw['ratio'][1] * sim) \
                    / sum(kw['ratio'])

                try:
                    iA = char_dict[charA]
                    iB = char_dict[charB]

                    # use the vowel scale
                    if charA[4] in self.vowels and charB[4] in self.vowels:
                        matrix[iA][iB] = matrix[iB][iA] = kw['vscale'] * rscore
                    else:
                        matrix[iA][iB] = matrix[iB][iA] = rscore
                except:
                    pass

        self.cscorer = misc.ScoreDict(self.chars, matrix)
        self._meta['scorer']['cscorer'] = self.cscorer

    def _get_partial_matrices(
            self,
            concept=False,
            method='sca',
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            mode='global',
            gop=-2,
            restriction='',
            **keywords
            ):
        """
        Function creates matrices for the purpose of partial cognate detection.
        """

        # set the defaults
        kw = dict(
            defaults=False,
            external_scorer=False,  # external scoring function
            imap_mode= False,
            sep=lingpy.settings.rcParams['morpheme_separator'],
            word_sep=lingpy.settings.rcParams['word_separator'],
            word_seps=lingpy.settings.rcParams['word_separators'],
            seps=lingpy.settings.rcParams['morpheme_separators'],
            tones='T',
            split_on_tones=False
        )
        kw.update(keywords)
        
        def function(idxA, idxB, sA, sB, **keywords):
            if method == 'lexstat':
                args = [
                        self[idxA, self._numbers][sA[0]:sA[1]],
                        self[idxB, self._numbers][sB[0]:sB[1]],
                        [self.cscorer[_charstring(
                            self[idxB, self._langid]
                            ), n]
                            for n in self[idxA, self._numbers][sA[0]:sA[1]]],
                        [self.cscorer[_charstring(
                            self[idxA, self._langid]), n]
                            for n in self[idxB, self._numbers][sB[0]:sB[1]]],
                        self[idxA, self._prostrings][sA[0]:sA[1]],
                        self[idxB, self._prostrings][sB[0]:sB[1]],
                        1,
                        scale,
                        factor,
                        self.cscorer,
                        mode,
                        restricted_chars,
                        1]
            elif method == 'sca':
                args = [
                        [n.split('.', 1)[1] for n in self[idxA,
                            self._numbers][sA[0]:sA[1]]],
                        [n.split('.', 1)[1] for n in self[idxB,
                            self._numbers][sB[0]:sB[1]]],
                        self[idxA, self._weights][sA[0]:sA[1]],
                        self[idxB, self._weights][sB[0]:sB[1]],
                        self[idxA, self._prostrings][sA[0]:sA[1]],
                        self[idxB, self._prostrings][sB[0]:sB[1]],
                        gop,
                        scale,
                        factor,
                        self.rscorer,
                        mode,
                        restricted_chars,
                        1]
            return calign.align_pair(*args)[2]
        
        concepts = [concept] if concept else sorted(self.rows)
        
        # we have two basic constraints in the algorithm:
        # a) set cognacy between morphemes in the same word to zero
        # b) set cognacy for those parts to zero which are superceded by
        # another part in all comparisons of two words
        # essentially, setting things to zero, means setting them to 1, since
        # we are dealing with distances here
        for c in concepts:
            
            indices = self.get_list(row=c, flat=True)
            matrix = []
            tracer = []
            
            # first assemble all partial parts
            trace = defaultdict(list) # stores where the stuff is in the matrix
            count = 0
            for idx in indices:
                
                # we need the slices for both words, so let's just take the
                # tokens for this time
                tokens = self[idx, self._segments]
                
                # now get the slices with the function
                slices = _get_slices(tokens, **kw)

                for i,slc in enumerate(slices):
                    tracer += [(idx, i, slc)]
                    trace[idx] += [(i, slc, count)]
                    count += 1
            
            if kw['imap_mode']:
                # now, iterate for each string pair, asses the scores, and make
                # sure, we only assign the best of those to the matrix

                matrix = [[0 for i in tracer] for j in tracer]
                # reset the self-constraints (we missed it before)


                for idxA, idxB in combinations(indices, r=2):
                    # iterate over all parts
                    scores = []
                    idxs = []
                    for i,sliceA,posA in trace[idxA]:
                        for j,sliceB,posB in trace[idxB]:
                            d = function(idxA, idxB, sliceA, sliceB)
                            scores += [d]
                            idxs += [(posA,posB)]
                    
                    visited_seqs = set([])
                    while scores:
                        min_score_index = scores.index(min(scores))
                        min_score = scores.pop(min_score_index)
                        posA, posB = idxs.pop(min_score_index)
                        if posA in visited_seqs or posB in visited_seqs:
                            matrix[posA][posB] = 1
                            matrix[posB][posA] = 1
                        else:
                            matrix[posA][posB] = min_score
                            matrix[posB][posA] = min_score
                            visited_seqs.add(posA)
                            visited_seqs.add(posB)
                for idx in indices:
                    for i,(_,sliceA,posA) in enumerate(trace[idx]):
                        for j,(_,sliceB,posB) in enumerate(trace[idx]):

                            if i < j:
                                matrix[posA][posB] = 1
                                matrix[posB][posA] = 1
            else:
                matrix = []
                for (idxA, posA, sliceA), (idxB, posB, sliceB) in combinations(tracer, r=2):
                    
                    if idxA == idxB:
                        d = 1
                    else:
                        try:
                            d = function(idxA, idxB, sliceA, sliceB)
                        except ZeroDivisionError:
                            lingpy.log.warning(
                                "Encountered Zero-Division for the comparison of "
                                "{0} and {1}".format(
                                    ''.join(self[idxA, self._tokens]),
                                    ''.join(self[idxB, self._tokens])))
                            d = 100
                    matrix += [d]
                matrix = lingpy.algorithm.misc.squareform(matrix)
            if not concept:
                yield c, tracer, matrix
            else:
                yield matrix

    def partial_cluster(
            self,
            method='sca',
            threshold=0.45,
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            mode='overlap',
            cluster_method='infomap',
            gop=-1,
            restriction='',
            ref='',
            external_function=None,
            split_on_tones=False,
            **keywords):
        """
        Cluster the words into partial cognate sets.

        Function for flat clustering of words into cognate sets.

        Parameters
        ----------
        method : {'sca','lexstat','edit-dist','turchin'} (default='sca')
            Select the method that shall be used for the calculation.
        cluster_method : {'upgma','single','complete', 'mcl'} (default='upgma')
            Select the cluster method. 'upgma' (:evobib:`Sokal1958`) refers to
            average linkage clustering, 'mcl' refers to the "Markov Clustering
            Algorithm" (:evobib:`Dongen2000`).
        threshold : float (default=0.3)
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
        inflation : {int, float} (default=2)
            Specify the inflation parameter for the use of the MCL algorithm.
        expansion : int (default=2)
            Specify the expansion parameter for the use of the MCL algorithm.
        
        """
        kw = dict(
                imap_mode=True,
                post_processing=True,
                inflation=2,
                expansion=2,
                max_steps=1000,
                add_self_loops=True,
                sep=lingpy.settings.rcParams['morpheme_separator'],
                word_sep=lingpy.settings.rcParams['word_separator'],
                word_seps=lingpy.settings.rcParams['word_separators'],
                seps=lingpy.settings.rcParams['morpheme_separators'],
                mcl_logs=lambda x: -np.log2((1 - x) ** 2)
                )
        kw.update(keywords)        

        # check for parameters and add clustering, in order to make sure that
        # analyses are not repeated
        if not hasattr(self, 'params'):
            self.params = {}
        self.params['partial_cluster'] = "{0}_{1}_{2:.2f}".format(
            method, cluster_method, threshold)
        self._stamp += '# Partial Cluster: ' + self.params['partial_cluster']

        matrices = self._get_partial_matrices(method=method, scale=scale,
                factor=factor, restricted_chars=restricted_chars, mode=mode,
                gop=gop, imap_mode=kw['imap_mode'],
                split_on_tones=split_on_tones)
        k = 0
        C = defaultdict(list) # stores the pcogids
        G = {} # stores the graphs
        with util.pb(desc='PARTIAL SEQUENCE CLUSTERING', total=len(self.rows)) as progress:
            for concept, trace, matrix in matrices:
                progress.update(1)
                lingpy.log.info('Analyzing concept {0}...'.format(concept))
                if external_function:
                    c = external_function(threshold, matrix,
                            taxa=list(range(len(matrix))), revert=True)
                elif cluster_method == 'infomap':
                    c = extra.infomap_clustering(threshold,
                            matrix, taxa=list(range(len(matrix))), 
                            revert=True)
                elif cluster_method == 'mcl':
                    c = clustering.mcl(threshold, matrix, 
                            taxa = list(range(len(matrix))),
                            max_steps=kw['max_steps'],
                            inflation=kw['inflation'],
                            expansion=kw['expansion'],
                            add_self_loops=kw['add_self_loops'],
                            logs=kw['mcl_logs'],
                            revert=True)
                elif cluster_method in ['upgma', 'single', 'complete', 'ward']:
                    c = clustering.flat_cluster(cluster_method,
                            threshold, matrix,
                            revert=True)
                else:
                    raise ValueError("No suitable cluster method specified.")
                
                for i, (idx, pos, slc) in enumerate(trace):
                    C[idx] += [c[i] + k]
                if kw['post_processing']:
                    _g = nx.Graph()
                    for i, (idx, pos, slc) in enumerate(trace):
                        _g.add_node((i,idx,pos))
                    remove_edges = []
                    for (i, n1), (j, n2) in util.combinations2(enumerate(_g.nodes())):
                        if C[n1[1]][n1[2]] == C[n2[1]][n2[2]]:
                            _g.add_edge(n1, n2)
                            if n1[1] == n2[1]:
                                # get scores for n1 and n2 with all the rest in
                                # the matrix to decide for one
                                sn1, sn2 = 0, 0
                                for i,row in enumerate(matrix):
                                    sn1 += matrix[i][n1[0]]
                                    sn2 += matrix[i][n2[0]]
                                sn1 = sn1 / len(matrix)
                                sn2 = sn2 / len(matrix)
                                if sn1 <= sn2:
                                    remove_edges += [n2]
                                else:
                                    remove_edges += [n1]
                    for node in remove_edges:
                        for edge in sorted(_g[node]):
                            _g.remove_edge(node, edge)

                    for i, coms in enumerate(nx.connected_components(_g)):
                        cogid = i + 1 + k
                        for j, idx, pos in coms:
                            C[idx][pos] = cogid
                    
                    G[concept] = _g

                k += len(matrix) + 1 
        self.add_entries(ref or self._partials, C, lambda x: x)
        self.graphs = G

    def add_cognate_ids(self, source, target, idtype='strict', override=False):
        """
        Compute normal cognate identifiers from partial cognate sets.

        Parameters
        ----------
        source: str
            Name of the source column in your wordlist file.
        target : str
            Name of the target column in your wordlist file.
        idtype : str (default="strict")
            Select between "strict" and "loose".
        override: bool (default=False)
            Specify whether you want to override existing columns.
        
        Notes
        -----
        While the computation of strict cognate IDs from partial cognate IDs is
        straightforward and just judges those words as cognate which are
        identical in all their parts, the computation of loose cognate IDs
        constructs a network between all words, draws lines between all words
        that share a common morpheme, and judges all connected components in this
        network as cognate.
        """
        if idtype == 'strict':
            
            tmp = defaultdict(list)
            for k in self._data:
                tmp[tuple(self[k, source])] += [k]
            idx = 1
            D = {}
            for vals in tmp.values():
                for k in vals:
                    D[k] = idx
                idx += 1
            self.add_entries(target, D, lambda x: x, override=override)
        elif idtype == 'loose':

            D = {}
            idx = 1
            for c in self.rows:
                idxs = self.get_list(row=c, flat=True)
                srcs = [self[k, source] for k in idxs]

                # get connected components
                g = nx.Graph()
                g.add_nodes_from(idxs)
                for (i, cogsA), (j, cogsB) in util.combinations2(zip(idxs, srcs)):
                     if [x for x in cogsA if x in cogsB]:
                         g.add_edge(i, j)
                for i,comps in enumerate(nx.connected_components(g)):
                    for comp in comps:
                        D[comp] = idx + i
                idx += (i+1)
            self.add_entries(target, D, lambda x: x, override=override)
        else:
            raise ValueError("The value you selected is not available.")
