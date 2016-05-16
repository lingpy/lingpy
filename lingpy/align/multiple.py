# *-* coding: utf-8 *-*
"""
Module provides classes and functions for multiple alignment analyses.
"""
from __future__ import print_function, division, unicode_literals
import logging
from itertools import combinations, combinations_with_replacement, product
from collections import defaultdict
from functools import partial

from lingpy.algorithm import calign
from lingpy.algorithm import talign
from lingpy.algorithm import cluster
from lingpy.algorithm import misc

from lingpy.thirdparty.cogent import LoadTree
from lingpy.sequence.sound_classes import (
    ipa2tokens, tokens2class, prosodic_string, prosodic_weights, pid,
)
from lingpy.settings import rcParams
from lingpy import log
from lingpy.util import setdefaults, identity, dotjoin, as_string 


class Multiple(object):
    """
    Basic class for multiple sequence alignment analyses.

    Parameters
    ----------
    seqs : list
        List of sequences that shall be aligned.

    Notes
    -----
    Depending on the structure of the sequences, further keywords can be
    specified that manage how the items get tokenized.

    """
    def __init__(self, seqs, **keywords):
        self.log = log.get_logger()
        # store input sequences, check whether tokens or strings are passed
        if isinstance(seqs[0], (list, tuple)):
            self.seqs = [' '.join(s) for s in seqs]
            self.tokens = [s for s in seqs]
        else:
            self.seqs = seqs
            self.tokens = []

        # define a tokenizer function for convenience
        kw = {
            "diacritics": rcParams['diacritics'],
            "vowels": rcParams['vowels'],
            "tones": rcParams['tones'],
            "combiners": rcParams['combiners'],
            "breaks": rcParams['breaks'],
            "stress": rcParams["stress"],
            "merge_vowels": rcParams["merge_vowels"],
            "unique_seqs": rcParams["unique_sequences"]
        }
        kw.update(keywords)

        self.numbers = []
        if self.tokens:
            for i, tokens in enumerate(self.tokens):
                self.numbers.append([dotjoin(i + 1, j + 1) for j in range(len(tokens))])
        else:
            # create a numerical representation of all sequences which reflects the
            # order of both their position and the position of their tokens. Before
            # this can be done, a tokenized version of all sequences has to be
            # created
            for i, seq in enumerate(self.seqs):
                # check for pre-tokenized strings
                tokens = ipa2tokens(seq, **kw)
                self.tokens.append(tokens)
                self.numbers.append([dotjoin(i + 1, j + 1) for j in range(len(tokens))])

        self.uniseqs = defaultdict(list)
        self.unique_seqs = kw["unique_seqs"]
        if self.unique_seqs:
            # create dictionary of all unique sequences, this is important, since
            # identical sequences should only be counted once in an alignment,
            # since they otherwise may disturb the analysis or slow it down
            for i, seq in enumerate(self.seqs):
                self.uniseqs[seq].append(i)
        else:
            # no uniqueness filtering
            self.uniseqs = range(0, len(self.seqs))

        self._length = len(self.uniseqs)

    def __len__(self):
        # the length of an alignment is defined as the number of unique
        # sequences present in the alignment
        return self._length

    def __str__(self):
        # if alignments are present, print the alignments
        # else, return all sequences
        lines = self.alm_matrix if self.alm_matrix else self.tokens
        return '\n'.join(['\t'.join(line) for line in lines])

    def __eq__(self, other):
        try:
            return self.alm_matrix == other.alm_matrix
        except:
            return False

    def __getitem__(self, idx):
        """
        Return specified values.
        """
        if isinstance(idx, tuple):
            if isinstance(idx[0], slice):
                return [x[idx[1]] for x in self.alm_matrix[idx[0]]]
            try:
                return self.alm_matrix[idx[0]][idx[1]]
            except:
                if idx[1] == 'w':
                    return self.seqs[idx[0]]
                if idx[1] == 'c':
                    return self.classes[idx[0]]
                if idx[1] == 't':
                    return self.tokens[idx[0]]
                if idx[1] == 'a':
                    return self.alm_matrix[idx[0]]
                return self.alm_matrix
        return self.seqs[idx]

    def _get(self, number, value='tokens', error=('X', '-')):
        """
        Method returns specific values of the class, depending on the index
        which is used.
        """
        # XXX this should be evaluated, maybe it is not needed in the future.
        if number == error[0]:
            return error[1]
        if number == '+':
            return "+"
        try:
            idxA, idxB = [int(i) - 1 for i in number.split('.')]

            if value == 'tokens':
                return self.tokens[idxA][idxB]
            if value == 'numbers':
                return self.numbers[idxA][idxB]
            if value == 'classes':
                return self.classes[idxA][idxB]
            if value == '_classes':
                return self._classes[idxA][idxB]
            if value == '_sonars':
                return self._sonars[idxA][idxB]
            if value == '_numbers':
                return self._numbers[idxA][idxB]
            if value == '_prosodics':
                return self._prosodics[idxA][idxB]
        except ValueError:
            if value == 'tokens':
                return self.tokens[int(number) - 1]
            if value == 'sonars':
                return self.sonars[int(number) - 1]
            if value == 'numbers':
                return self.numbers[int(number) - 1]
            if value == 'classes':
                return self.classes[int(number) - 1]
            if value == '_sonars':
                return self._sonars[int(number) - 1]
            if value == '_numbers':
                return self._numbers[int(number) - 1]
            if value == '_classes':
                return self._classes[int(number) - 1]

    def _set_model(
            self, model=None, classes=True, sonar=True, sonars=False, scoredict={}):
        """
        Method defines a specific class model for the calculation.

        Parameters
        ----------
        model : { None ~lingpy.data.model.Model } (default=None)
            A sound class model.
        """
        # check whether model is a string
        if not hasattr(model, 'name'):
            model = rcParams[model]

        # check for keyword classes
        if not classes:
            classify = identity
        else:
            self.model = model or rcParams['sca']
            classify = lambda x: tokens2class(x, self.model)

        # create the sound-classes or the fake classes
        self.classes = [cls for cls in map(classify, self.tokens)]

        # once a class model is defined, there may be identical sequences,
        # which in IPA terms are different. In order to avoid computing
        # alignments for these identical sequences, a dictionary is created
        # which stores references to all identical sequences, thus allowing to
        # compute only one alignment for each set of identical sequences
        indices = defaultdict(list)
        for i, seq in enumerate(self.classes):
            indices[tuple(seq)].append(i)

        # create additional matrices for the internal representation of the
        # class sequences
        if self.unique_seqs:
            keys = [val[0] for val in indices.values()]
        else:
            keys = range(len(self.classes))
        self.height = len(keys)

        # add the classes
        self._classes = [self.classes[key] for key in keys]
        self._numbers = [[dotjoin(i + 1, j + 1) for j in
                          range(len(self._classes[i]))] for i in range(self.height)]

        # create an index which allows to quickly interchange between classes
        # and given sequences (trivial without sequence uniqueness
        if self.unique_seqs:
            self.int2ext = {i: indices[tuple(self._classes[i])] for i in range(len(keys))}
        else:
            self.int2ext = {i: [i] for i in range(len(keys))}

        # -> # create external to internal in order to allow for a quick switching
        # -> # of the vals
        # -> self.ext2int = {}
        # -> for k,vals in self.int2ext.items():
        # ->     for v in vals:
        # ->         self.ext2int[v] = k

        # store sonars if they are passed as a list
        if sonar and sonars:  # == list:
            self._sonars = [sonars[key] for key in keys]
            # -> self._sonars = [0 for i in range(len(sonar))]
            # -> for i in range(len(self._sonars)):
            # ->     self._sonars[self.ext2int[i]] = sonar[i]
            self._prostrings = list([prosodic_string(s) for s in self._sonars])
        # create sonars if the argument is true
        elif sonar:
            self._sonars = list(
                map(lambda x: [int(t) for t in tokens2class(
                    # XXX change this part
                    x, rcParams['art'], stress=rcParams['stress'])],
                    [self.tokens[key] for key in keys]))
            if log.get_level() <= logging.DEBUG:
                for _i, _sonar in enumerate(self._sonars):
                    if 0 in _sonar:
                        self.log.warn(
                            "Sequence {0} contains unrecognized characters!".format(
                                self.seqs[self.int2ext[_i][0]]))
            self._prostrings = list([prosodic_string(s) for s in self._sonars])
        # do nothing if no arguments are passed
        else:
            self._sonars = False
            self._prostrings = False

        # create a scoredict for the calculation of alignment analyses
        # append the scorer if it is given with the model
        def scorer(x, y):
            if classes:
                return self.model.scorer[x, y]
            if scoredict:
                return scoredict[x, y]
            return 1.0 if x == y else -1.0

        self.scoredict = {}
        for (i, seqA), (j, seqB) in combinations_with_replacement(
                enumerate(self._numbers), 2):
            if i < j:
                for (numA, numB) in product(seqA, seqB):
                    self.scoredict[numA, numB] = scorer(
                        self._get(numA, '_classes'), self._get(numB, '_classes'))
                    self.scoredict[numB, numA] = self.scoredict[numA, numB]
            elif i == j:
                for num in seqA:
                    char = self._get(num, '_classes')
                    self.scoredict[num, num] = scorer(char, char)

    def _set_scorer(self, score_mode='classes'):
        """
        Functions sets the scorer to the simple class model or to the library
        model.
        """
        if score_mode == 'classes':
            self.scorer = self.scoredict
        elif score_mode == 'library':
            self.scorer = self.library

    def _get_pairwise_alignments(
            self,
            mode='global',
            gop=-2,
            scale=0.5,
            factor=0.3,
            restricted_chars='T_',
            **keywords):
        """
        Function calculates all pairwise alignments from the data.
        """
        if 'transform' not in keywords:
            keywords['transform'] = rcParams['align_transform']

        # create array for alignments
        self._alignments = [[0 for i in range(self.height)] for i in range(self.height)]

        # create the distance matrix
        self.matrix = []

        # check for the mode, if sonority profiles are not chose, take the
        # simple alignment function
        if self._sonars:
            make_pro_weights = partial(prosodic_weights, _transform=keywords['transform'])

            # get the weights
            if not hasattr(self, 'weights'):
                self._weights = list(map(make_pro_weights, self._prostrings))

            alignments = calign.align_pairwise(
                self._numbers,
                self._weights,
                self._prostrings,
                gop,
                scale,
                factor,
                self.scorer,
                restricted_chars,
                mode)
            k = 0
            for i, j in combinations_with_replacement(range(self.height), 2):
                almA, almB, sim, dist = alignments[k]
                k += 1
                if i < j:
                    if mode == 'local':
                        almA = almA[1]
                        almB = almB[1]
                    self._alignments[i][j] = [almA, almB, sim]
                    self._alignments[j][i] = [almA, almB, sim]
                    self.matrix += [dist]
                elif i == j:
                    self._alignments[i][j] = [almA, almB, sim]
        else:
            alignments = talign.align_pairwise(
                self._numbers, gop, scale, self.scorer, mode)
            k = 0
            for i, j in combinations_with_replacement(range(self.height), 2):
                almA, almB, sim, dist = alignments[k]
                k += 1
                if i < j:
                    if mode == 'local':
                        almA = almA[1]
                        almB = almB[1]
                    self._alignments[i][j] = [almA, almB, sim]
                    self._alignments[j][i] = [almA, almB, sim]
                    self.matrix += [dist]
                elif i == j:
                    self._alignments[i][j] = [almA, almB, sim]

        self.matrix = misc.squareform(self.matrix)

    def _create_library(self):
        """
        Method creates an extended library for alignments using the Tcoffee
        approach.
        """
        self.library = {}

        # create library for non-sound-class approaches
        if not self._sonars:
            for numA, numB in combinations_with_replacement(self._numbers, 2):
                for k, l in product(numA, numB):
                    self.library[k, l] = 0.0
                    self.library[l, k] = 0.0
        else:
            # note that we somehow HAVE to include a sensitivity for V-C
            # distinctions in the library mode, otherwise it may get complicated
            # sometimes, therefore, the library is initialized by setting only the
            # scores for c-c and v-v matches to 0, the other scores get their
            # original penalty defined by the old scorer
            for (i, numA), (j, numB) in combinations_with_replacement(
                    enumerate(self._numbers), 2):
                if i < j:
                    for k, l in product(numA, numB):
                        # see the comment above for the add-on in this
                        # line
                        a = self._get(k, '_sonars')
                        b = self._get(l, '_sonars')
                        if a >= 7 or b >= 7 and a + b < 14:
                            self.library[k, l] = self.scoredict[k, l]
                            self.library[l, k] = self.scoredict[l, k]
                        else:
                            self.library[k, l] = 0.0
                            self.library[l, k] = 0.0
                elif i == j:
                    for k, l in product(numA, numB):
                        self.library[k, l] = 0.0
                        self.library[l, k] = 0.0

    def _extend_library(self):
        """
        Extend the library by new alignments.
        """
        # add the residue-pairs of all aligned sequences first
        for i, j in [(i, j) for i in range(self.height)
                     for j in range(self.height) if i <= j]:
            for m, n in zip(self._alignments[i][j][0], self._alignments[i][j][1]):
                if m != "-" and n != "-":
                    # add the values to the library
                    # the similarity score is determined by adding taking the
                    # average of matrix score and the similarity score of the
                    # alignment of both sequences
                    score = self.scorer[m, n]
                    sim = self._alignments[i][j][2] / float(len(self._alignments[i][j][0]))
                    self.library[m, n] += (sim + score) / 2.0
                    self.library[n, m] = self.library[m, n]

        # add the residue-pairs resulting from an alignment via a third sequence

        # create the indices for the loop
        mappings = (
            (i, j, k)
            for i in range(self.height)
            for j in range(self.height)
            for k in range(self.height)
            if i <= j and k != i and k != j)

        for i, j, k in mappings:
            almI, almIK, simIK = self._alignments[i][k]
            almJ, almJK, simJK = self._alignments[j][k]

            # determine, which of the values occur in both alignments
            # with the third sequence
            for char in self._numbers[k]:
                try:
                    valI = almI[almIK.index(char)]
                    valJ = almJ[almJK.index(char)]
                    if valI != "-" and valJ != "-":
                        score = self.scorer[valI, valJ]
                        sim = min(simIK, simJK) / ((len(almIK) + len(almJK)) / 2.0)

                        self.library[valI, valJ] += (sim + score) / 2.0
                        self.library[valJ, valI] = self.library[valI, valJ]
                except:
                    pass

    def _make_guide_tree(self, tree_calc='upgma'):
        """
        Create the guide tree using either the UPGMA or the Neighbor-Joining
        algorithm.
        """
        clusters = {i[0]: [i[1]] for i in zip(range(self.height), range(self.height))}

        # create the tree matrix
        self.tree_matrix = []

        # carry out the clustering
        if tree_calc == 'upgma':
            cluster._upgma(clusters, self.matrix, self.tree_matrix)
        elif tree_calc == 'neighbor':
            cluster._neighbor(clusters, self.matrix, self.tree_matrix)
        else:
            raise ValueError(
                'Method <' + tree_calc + '> for tree calculation not available.')

        # create a newick-representation of the string
        self.tree = LoadTree(cluster._tree2nwk(
            self.tree_matrix, [''.join(str(c)) for c in self._classes], False))

    def _align_profile(
            self,
            almsA,
            almsB,
            mode='global',
            gop=-3,
            scale=0.5,
            factor=0,
            gap_weight=0.5,
            return_similarity=False,
            iterate=False,
            restricted_chars="T_"):
        profileA = misc.transpose(almsA)
        profileB = misc.transpose(almsB)

        # calculate profile length and profile depth for both profiles
        o = len(profileA[0])
        p = len(profileB[0])

        # create the weights by which the gap opening penalties will be modified
        sonarA = [[self._get(char, value='_sonars', error=('X', 0))
                   for char in line] for line in profileA]
        sonarB = [[self._get(char, value='_sonars', error=('X', 0))
                   for char in line] for line in profileB]

        # get the consensus string for the sonority profiles
        try:
            consA = [
                int(sum([k for k in col if k != 0]) /
                    len([k for k in col if k != 0]) + 0.5) for col in sonarA]
            consB = [
                int(sum([k for k in col if k != 0]) /
                    len([k for k in col if k != 0]) + 0.5) for col in sonarB]
        except:
            try:
                consA = [
                    int(sum([k for k in col if k >= 0]) /
                        len([k for k in col if k >= 0]) + 0.5) for col in sonarA]
                consB = [
                    int(sum([k for k in col if k >= 0]) /
                        len([k for k in col if k >= 0]) + 0.5) for col in sonarB]
                self.log.warn("There are empty segments in the consensus.")
                self.log.info(
                    '',
                    extra=dict(lines=[' '.join([str(x) for x in cons])
                                      for cons in [consA, consB]]))
            except:
                self.log.error(
                    "Failed to compute the consensus string.",
                    extra=dict(lines=[
                        sonarA, sonarB,
                        almsA[0], [self._get(n_, 'tokens') for n_ in almsA[0]],
                        almsB[0], [self._get(n_, 'tokens') for n_ in almsB[0]]
                    ]))

        prosA = prosodic_string(consA)
        prosB = prosodic_string(consB)

        self.log.debug('', extra=dict(lines=[(prosA, consA), (prosB, consB)]))
        weightsA, weightsB = prosodic_weights(prosA), prosodic_weights(prosB)

        # carry out the alignment
        almA, almB, sim = calign.align_profile(
            profileA,
            profileB,
            weightsA,
            weightsB,
            prosA,
            prosB,
            gop,
            scale,
            factor,
            self.scorer,
            restricted_chars,
            mode,
            gap_weight)

        if return_similarity:
            return sim

        # trace the gaps inserted in both aligned profiles and insert them
        # in the original profiles
        for i in range(len(almA)):
            if almA[i] == '-':
                profileA.insert(i, o * ['X'])
            elif almB[i] == '-':
                profileB.insert(i, p * ['X'])

        # invert the profiles and the weight matrices by turning columns
        # into rows and rows into columns
        profileA = misc.transpose(profileA)
        profileB = misc.transpose(profileB)

        # return the aligned profiles and weight matrices
        if iterate:
            return profileA, profileB

        return profileA + profileB

    def _talign_profile(
            self,
            almsA,
            almsB,
            mode='global',
            gop=-3,
            scale=0.5,
            gap_weight=0.5,
            return_similarity=False,
            iterate=False):
        """
        Align profiles for tokens, not sound classes.
        """
        profileA = misc.transpose(almsA)
        profileB = misc.transpose(almsB)

        # calculate profile length and profile depth for both profiles
        o = len(profileA[0])
        p = len(profileB[0])

        # carry out the alignment
        almA, almB, sim = talign.align_profile(
            profileA, profileB, gop, scale, self.scorer, mode, gap_weight)

        if return_similarity:
            return sim

        # trace the gaps inserted in both aligned profiles and insert them
        # in the original profiles
        for i in range(len(almA)):
            if almA[i] == '-':
                profileA.insert(i, o * ['X'])
            elif almB[i] == '-':
                profileB.insert(i, p * ['X'])

        # invert the profiles and the weight matrices by turning columns
        # into rows and rows into columns
        profileA = misc.transpose(profileA)
        profileB = misc.transpose(profileB)

        # return the aligned profiles and weight matrices
        if iterate:
            return profileA, profileB

        return profileA + profileB

    def _merge_alignments(
            self,
            mode='global',
            gop=-3,
            scale=0.5,
            factor=0,
            gap_weight=0.5,
            restricted_chars='T_'):
        # create the lists which will store the current stages of the
        # alignment process
        seq_ord = [[i] for i in range(self.height)]
        alm_lst = [[seq] for seq in self._numbers[:]]

        # start the iteration through the tree array: the first two lines
        # in the matrix contain the ids of the sequences in the array,
        # which are aligned along the tree
        if self._sonars:
            algorithm = self._align_profile
            kw = dict(factor=factor, restricted_chars=restricted_chars)
        else:
            algorithm = self._talign_profile
            kw = {}

        for row in self.tree_matrix:
            m, n = int(row[0]), int(row[1])
            seq_ord.append(seq_ord[m] + seq_ord[n])

            alms = algorithm(
                alm_lst[m],
                alm_lst[n],
                mode=mode,
                gop=gop,
                scale=scale,
                gap_weight=gap_weight,
                **kw)
            alm_lst.append(alms)

        # get the last stage of each alignment process
        alm_lst = alm_lst[-1]

        # restore the original order of the strings in the alignment
        sorter = seq_ord[-1][:]
        sorter.reverse()
        alm_lst = sorted(alm_lst, key=lambda x: sorter.pop())

        # create the matrix which stores all alignments
        self._alm_matrix = alm_lst

        # calculate the sonority profile
        if self._sonars:
            tmp = misc.transpose(alm_lst)
            sonars = [[self._get(char, value='_sonars', error=('X', 0))
                       for char in line] for line in tmp]
            try:
                consensus = [
                    int(sum([k for k in col if k != 0]) /
                        len([k for k in col if k != 0]) + 0.5) for col in sonars]
            except:
                try:
                    consensus = [
                        int(sum([k for k in col if k >= 0]) /
                            len([k for k in col if k >= 0]) + 0.5) for col in sonars]
                    self.log.warn("There are empty segments in the consensus.")
                    self.log.info('', extra=dict(lines=[consensus]))
                except:
                    consensus = []
                    self.log.error("Failed to compute the consensus string.")
            self._sonority_consensus = consensus

    def _update_alignments(self):
        self.alm_matrix = [0 for i in range(len(self.numbers))]

        for i, line in enumerate(self._alm_matrix):
            indices = self.int2ext[i]
            for j in indices:
                numbers = []
                for num in line:
                    try:
                        numbers.append(dotjoin(j + 1, num.split('.')[1]))
                    except:
                        numbers.append('X')
                self.alm_matrix[j] = [self._get(num, 'tokens') for num in numbers]

    def prog_align(self, **keywords):
        """
        Carry out a progressive alignment analysis of the input sequences.

        Parameters
        ----------

        model : { "dolgo", "sca", "asjp" } (defaul="sca")
            A string indicating the name of the :py:class:`Model \
            <lingpy.data.model>` object that shall be used for the analysis.
            Currently, three models are supported:

            * "dolgo" -- a sound-class model based on :evobib:`Dolgopolsky1986`,

            * "sca" -- an extension of the "dolgo" sound-class model based on
              :evobib:`List2012b`, and

            * "asjp" -- an independent sound-class model which is based on the
              sound-class model of :evobib:`Brown2008` and the empirical data
              of :evobib:`Brown2011` (see the description in
              :evobib:`List2012`.

        mode : { "global", "dialign" } (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

        gop : int (default=-2)
            The gap opening penalty (GOP) used in the analysis.

        scale : float (default=0.5)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1982`.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.

        tree_calc : { "neighbor", "upgma" } (default="upgma")
            The cluster algorithm which shall be used for the calculation of
            the guide tree. Select between ``neighbor``, the Neighbor-Joining
            algorithm (:evobib:`Saitou1987`), and ``upgma``, the UPGMA
            algorithm (:evobib:`Sokal1958`).

        guide_tree : tree_matrix
            Use a custom guide tree instead of performing a cluster algorithm
            for constructing one based on the input similarities. The use of this
            option makes the tree_calc option irrelevant.

        gap_weight : float (default=0.5)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        restricted_chars : string (default="T")
            Define which characters of the prosodic string of a sequence
            reflect its secondary structure (cf. :evobib:`List2012b`) and
            should therefore be aligned specifically. This defaults to "T",
            since this is the character that represents tones in the prosodic
            strings of sequences.

        """
        # set up the defaults parameters stored in the kw dictionary
        kw = dict(
            model=rcParams['model'],
            mode=rcParams['align_mode'],
            scale=rcParams['align_scale'],
            factor=rcParams['align_factor'],
            tree_calc=rcParams['align_tree_calc'],
            restricted_chars=rcParams['restricted_chars'],
            classes=rcParams['align_classes'],
            sonar=rcParams['align_sonar'],
            sonars=False,
            scoredict=rcParams['align_scorer'],
            gop=rcParams['align_gop'],
            gap_weight=rcParams['align_gap_weight']
        )
        kw.update(keywords)

        # fixing a but to avoid that defining models as string will yield an error
        if not hasattr(kw['model'], 'name'):
            kw['model'] = rcParams[kw['model']]

        # define the model for convenience
        model = kw['model']

        # create a string with the current parameters
        self.params = '_'.join([
            'prog',
            model.name,
            str(kw['gop']),
            '{0:.1f}'.format(kw['scale']),
            '{0:.1f}'.format(kw['factor']),
            kw['tree_calc'],
            '{0:.1f}'.format(kw['gap_weight']),
            kw['restricted_chars']
        ])

        self._set_model(model, kw['classes'], kw['sonar'], kw['sonars'], kw['scoredict'])
        self._set_scorer('classes')

        self._get_pairwise_alignments(
            gop=kw['gop'],
            scale=kw['scale'],
            factor=kw['factor'],
            restricted_chars=kw['restricted_chars'])

        if 'guide_tree' in kw.keys():
            self.tree_matrix = kw['guide_tree']
        else:
            self._make_guide_tree(tree_calc=kw['tree_calc'])

        self._merge_alignments(
            mode=kw['mode'],
            gop=kw['gop'],
            scale=kw['scale'],
            factor=kw['factor'],
            restricted_chars=kw['restricted_chars'],
            gap_weight=kw['gap_weight'])

        self._update_alignments()

    def lib_align(self, **keywords):
        """
        Carry out a library-based progressive alignment analysis of the sequences.

        Notes
        -----
        In contrast to traditional progressive multiple sequence alignment
        approaches such as :evobib:`Feng1981` and :evobib:`Thompson1994`,
        library-based progressive alignment :evobib:`Notredame2000` is based on
        a pre-processing of the data where the information given in global and
        local pairwise alignments of the input sequences is used to derive a
        refined scoring function (*library*) which is later used in the
        progressive phase.


        Parameters
        ----------

        model : { "dolgo", "sca", "asjp" } (default="sca")
            A string indicating the name of the :py:class:`Model \
            <lingpy.data.model>` object that shall be used for the analysis.
            Currently, three models are supported:

            * "dolgo" -- a sound-class model based on :evobib:`Dolgopolsky1986`,

            * "sca" -- an extension of the "dolgo" sound-class model based on
              :evobib:`List2012b`, and

            * "asjp" -- an independent sound-class model which is based on the
              sound-class model of :evobib:`Brown2008` and the empirical data
              of :evobib:`Brown2011` (see the description in
              :evobib:`List2012`.

        mode : { "global", "dialign" } (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

        modes : list (default=[("global",-10,0.6),("local",-1,0.6)])
            Indicate the mode, the gap opening penalties (GOP), and the gap extension
            scale (GEP scale), of the pairwise alignment analyses which
            are used to create the library.

        gop : int (default=-5)
            The gap opening penalty (GOP) used in the analysis.

        gep_scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1982`.

        factor : float (default=1)
            The factor by which the initial and the descending position shall
            be modified.

        tree_calc : { "neighbor", "upgma" } (default="upgma")
            The cluster algorithm which shall be used for the calculation of
            the guide tree. Select between ``neighbor``, the Neighbor-Joining
            algorithm (:evobib:`Saitou1987`), and ``upgma``, the UPGMA
            algorithm (:evobib:`Sokal1958`).

        guide_tree : tree_matrix
            Use a custom guide tree instead of performing a cluster algorithm
            for constructing one based on the input similarities. The use of this
            option makes the tree_calc option irrelevant.

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        restricted_chars : string (default="T")
            Define which characters of the prosodic string of a sequence
            reflect its secondary structure (cf. :evobib:`List2012b`) and
            should therefore be aligned specifically. This defaults to "T",
            since this is the character that represents tones in the prosodic
            strings of sequences.

        """
        # set up the defaults parameters stored in the kw dictionary
        kw = dict(
            model=rcParams['model'],
            mode=rcParams['align_mode'],
            modes=rcParams['align_modes'],
            scale=rcParams['align_scale'],
            factor=rcParams['align_factor'],
            tree_calc=rcParams['align_tree_calc'],
            restricted_chars=rcParams['restricted_chars'],
            classes=rcParams['align_classes'],
            sonar=rcParams['align_sonar'],
            scoredict=rcParams['align_scorer'],
            gop=rcParams['align_gop'],
            gap_weight=rcParams['align_gap_weight'],
            sonars=False)
        kw.update(keywords)

        # fixing a but to avoid that defining models as string will yield an error
        if not hasattr(kw['model'], 'name'):
            kw['model'] = rcParams[kw['model']]

        # create a string with the current parameters
        params = [
            'lib',
            kw['model'].name,
            kw['mode'],
            '{0:.1f}'.format(kw['factor']),
            kw['tree_calc'],
            '{0:.1f}'.format(kw['gap_weight']),
            kw['restricted_chars']
        ]

        # append parameters to the params-string
        mode_params = []
        for m in kw['modes']:
            mode_params.append('{0[0]}x{0[1]}x{0[2]:.2f}'.format(m))

        mode_params = 'y'.join(mode_params)

        params[3] = mode_params
        self.params = '_'.join(params)

        self._set_model(
            kw['model'], kw['classes'], kw['sonar'], kw['sonars'], kw['scoredict'])

        self._set_scorer('classes')

        # start to create the library, note that scales and factors are set to
        # zero here, since scales and zeros are only useful in
        # profile-alignments. they eventually disturb pairwise alignments,
        # which is why it is important to keep their influence low when
        # creating the library from pairwise alignments
        self._create_library()
        for run in kw['modes']:
            self._get_pairwise_alignments(
                run[0], run[1], run[2], kw['factor'], kw['restricted_chars'])
            self._extend_library()

        self._set_scorer('library')
        self._get_pairwise_alignments(
            kw['mode'], 0, 0.0, kw['factor'], kw['restricted_chars'])

        if 'guide_tree' in kw.keys():
            self.tree_matrix = kw['guide_tree']
        else:
            self._make_guide_tree(tree_calc=kw['tree_calc'])

        # merge the alignments, not that the scale doesn't really influence any
        # of the results here, since gap scores are set to 0, gapping should be
        # the same in all positions, the factor, however, eventually influences
        # the score, since it changes character mappings as well
        self._merge_alignments(
            kw['mode'], 0, 0.0, 0, kw['gap_weight'], kw['restricted_chars'])

        self._update_alignments()

    def _reduce_gap_sites(self, msa, gap='X'):
        """
        Method reduces all columns from an MSA when there are only gaps. This
        method is important for the iterative procedures.
        """
        # XXX new_msa = np.array(msa[:])
        new_msa = [m for m in msa]
        no_gap_index = []
        for i in range(len(new_msa[0])):
            # XXX if list(new_msa[:,i]).count(gap) != len(new_msa):
            if [line[i] for line in new_msa].count(gap) != len(new_msa):
                no_gap_index.append(i)

        new_msa = [[line[i] for i in no_gap_index] for line in new_msa]
        # XXX new_msa = new_msa[:,no_gap_index].tolist()
        return new_msa

    def _split(self, idx):
        """
        Split an MSA into two parts and retain their indices.
        """
        # XXX
        # create the inverted index
        idxA = idx
        idxB = [i for i in range(self.height) if i not in idx]

        # get idxA
        almA = [self._alm_matrix[i] for i in idxA]
        almB = [self._alm_matrix[i] for i in idxB]

        partA = self._reduce_gap_sites(almA)
        partB = self._reduce_gap_sites(almB)

        # XXX partA = self._reduce_gap_sites(self._alm_matrix[idxA])
        # XXX partB = self._reduce_gap_sites(self._alm_matrix[idxB])

        return partA, partB, idxA, idxB

    def _join(self, almA, almB, idxA, idxB):
        """
        Join two aligned MSA by their index.
        """
        m = len(almA[0])
        out_alm = [[0 for i in range(m)] for j in range(self.height)]

        for idx, alm in [(idxA, almA), (idxB, almB)]:
            for i in range(len(alm)):
                out_alm[idx[i]] = alm[i]

        # XXX out_alm = np.array(out_alm)
        return out_alm

    def _iter(
            self,
            idx_list,
            mode='global',
            gop=-3,
            scale=0.5,
            factor=0.0,
            gap_weight=0.5,
            check='final',
            restricted_chars='T_'):
        """
        Split an MSA into two parts and realign them.
        """
        sop = self.sum_of_pairs(gap_weight=gap_weight)
        alm_matrix = [[cell for cell in line] for line in self._alm_matrix]
        # XXX .copy()

        if len(idx_list) == 1:
            return

        for idx in idx_list:
            almA, almB, idxA, idxB = self._split(idx)
            almA, almB = self._align_profile(
                almA,
                almB,
                mode=mode,
                iterate=True,
                gop=gop,
                scale=scale,
                factor=factor,
                gap_weight=gap_weight)
            new_alm = self._join(almA, almB, idxA, idxB)

            self._alm_matrix = new_alm
            if check == 'immediate':
                new_sop = self.sum_of_pairs()
                if new_sop < sop:
                    self._alm_matrix = alm_matrix
                else:
                    sop = new_sop

        if check == 'final':
            new_sop = self.sum_of_pairs(gap_weight=gap_weight)
            if new_sop < sop:
                self._alm_matrix = alm_matrix

        self._update_alignments()

    def sum_of_pairs(self, alm_matrix='self', mat=None, gap_weight=0.0, gop=-1):
        """
        Calculate the sum-of-pairs score for a given alignment analysis.

        Parameters
        ----------

        alm_matrix : { "self", "other" } (default="self")
            Indicate for which MSA the sum-of-pairs score shall be calculated.

        mat : { None, list }
            If "other" is chosen as an option for **alm_matrix**, define for
            which matrix the sum-of-pairs score shall be calculated.

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        Returns
        -------
        The sum-of-pairs score of the alignment.
        """
        if alm_matrix == 'self':
            alm_matrix = self._alm_matrix
        else:
            alm_matrix = mat

        lenM = len(alm_matrix[0])

        score = 0.0
        args, kw = [], {}
        if self._sonars:
            algorithm = calign
            kw = dict(gap_weight=gap_weight)
        else:
            algorithm = talign
            args = [gop, gap_weight]

        for i in range(lenM):
            score += algorithm.score_profile(
                [line[i] for line in alm_matrix],
                [line[i] for line in alm_matrix],
                self.scorer,
                *args,
                **kw)
        return score / lenM

    def _swap_sum_of_pairs(self, alm_matrix, gap_weight=1.0, swap_penalty=-5):
        lenM = len(alm_matrix[0])
        score = 0.0
        algorithm = calign if self._sonars else talign

        for i in range(lenM):
            score += algorithm.swap_score_profile(
                [line[i] for line in alm_matrix],
                [line[i] for line in alm_matrix],
                self.scorer,
                gap_weight=gap_weight,
                swap_penalty=swap_penalty)

        return score / lenM

    def iterate_orphans(
            self,
            check='final',
            mode='global',
            gop=-3,
            scale=0.5,
            factor=0,
            gap_weight=1.0,
            restricted_chars='T_'):
        """
        Iterate over the most divergent sequences in the sample.

        Parameters
        ----------

        check : string (default="final")
            Specify when to check for improved sum-of-pairs scores: After each
            iteration ("immediate") or after all iterations have been carried
            out ("final").

        mode : { "global", "overlap", "dialign" } (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

            * "overlap" -- semi-global alignment, where gaps introduced in the
              beginning and the end of a sequence do not score.

        gop : int (default=-5)
            The gap opening penalty (GOP) used in the analysis.

        gep_scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1981`.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        Notes
        -----
        The most divergent sequences are those whose average distance to all
        other sequences is above the average distance of all sequence pairs.

        See also
        --------
        Multiple.iterate_clusters
        Multiple.iterate_similar_gap_sites
        Multiple.iterate_all_sequences

        """
        orphans = []
        means = [sum(line) / len(line) for line in self.matrix]  # XXX self.matrix.mean()
        means = sum(means) / len(means)

        for i, line in enumerate(self.matrix):
            if sum(line) / len(line) > means:
                orphans.append([i])

        self._iter(
            orphans,
            check=check,
            mode=mode,
            scale=scale,
            gop=gop,
            factor=factor,
            gap_weight=gap_weight,
            restricted_chars=restricted_chars)

    def iterate_clusters(
            self,
            threshold,
            check='final',
            mode='global',
            gop=-3,
            scale=0.5,
            factor=0,
            gap_weight=1,
            restricted_chars='T_'):
        """
        Iterative refinement based on a flat cluster analysis of the data.

        Notes
        -----
        This method uses the :py:func:`lingpy.algorithm.clustering.flat_upgma`
        function in order to retrieve a flat cluster of the data.

        Parameters
        ----------

        threshold : float
            The threshold for the flat cluster analysis.

        check : string (default="final")
            Specify when to check for improved sum-of-pairs scores: After each
            iteration ("immediate") or after all iterations have been carried
            out ("final").

        mode : { "global", "overlap", "dialign" } (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * 'global' -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * 'dialign' -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

            * 'overlap' -- semi-global alignment, where gaps introduced in the
              beginning and the end of a sequence do not score.

        gop : int (default=-5)
            The gap opening penalty (GOP) used in the analysis.

        gep_scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1981`.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        See also
        --------
        Multiple.iterate_similar_gap_sites
        Multiple.iterate_all_sequences

        """
        # don't calculate this if there are less than 5 sequences
        if len(self.seqs) < 3:
            return

        # create the clusters
        clusters = {i[0]: [i[1]] for i in zip(range(self.height), range(self.height))}

        cluster._flat_upgma(clusters, self.matrix, threshold)
        self._iter(
            clusters.values(),
            check=check,
            mode=mode,
            scale=scale,
            gop=gop,
            factor=factor,
            gap_weight=gap_weight,
            restricted_chars=restricted_chars)

    def iterate_similar_gap_sites(
            self,
            check='final',
            mode='global',
            gop=-3,
            scale=0.5,
            factor=0,
            gap_weight=1,
            restricted_chars='T_'):
        """
        Iterative refinement based on the *Similar Gap Sites* heuristic.

        Notes
        -----
        This heuristic is fairly simple. The idea is to try to split a given
        MSA into partitions with identical gap sites.

        Parameters
        ----------

        check : { "final", "immediate" } (default="final")
            Specify when to check for improved sum-of-pairs scores: After each
            iteration ("immediate") or after all iterations have been carried
            out ("final").

        mode : { "global", "overlap", "dialign" } (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * 'global' -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * 'dialign' -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

            * 'overlap' -- semi-global alignment, where gaps introduced in the
              beginning and the end of a sequence do not score.

        gop : int (default=-5)
            The gap opening penalty (GOP) used in the analysis.

        gep_scale : float (default=0.5)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1982`.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.

        gap_weight : float (default=1)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When, e.g., set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        See also
        --------
        Multiple.iterate_clusters
        Multiple.iterate_all_sequences
        Multiple.iterate_orphans

        """
        self._similar_gap_sites()

        if len(self.gap_dict) == 1:
            return
        self._iter(
            list(self.gap_dict.values()),
            check=check,
            mode=mode,
            scale=scale,
            gop=gop,
            factor=factor,
            gap_weight=gap_weight)

    def iterate_all_sequences(
            self,
            check="final",
            mode="global",
            gop=-3,
            scale=0.5,
            factor=0,
            gap_weight=1,
            restricted_chars="T_"):
        """
        Iterative refinement based on a complete realignment of all sequences.

        Notes
        -----
        This method essentially follows the iterative method of
        :evobib:`Barton1987` with the exception that an MSA has already been
        calculated.

        Parameters
        ----------

        check : { "final", "immediate" } (default="final")
            Specify when to check for improved sum-of-pairs scores: After each
            iteration ("immediate") or after all iterations have been carried
            out ("final").

        mode : { "global", "overlap", "dialign" } (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

            * "overlap" -- semi-global alignment, where gaps introduced in the
              beginning and the end of a sequence do not score.

        gop : int (default=-5)
            The gap opening penalty (GOP) used in the analysis.

        gep_scale : float (default=0.5)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1981`.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        See also
        --------
        Multiple.iterate_clusters
        Multiple.iterate_similar_gap_sites
        Multiple.iterate_orphans

        """
        self._iter(
            [[i] for i in range(self.height)],
            check=check,
            mode=mode,
            scale=scale,
            gop=gop,
            factor=factor,
            gap_weight=gap_weight,
            restricted_chars=restricted_chars)

    def get_peaks(self, gap_weight=0):
        """
        Calculate the profile score for each column of the alignment.

        Parameters
        ----------

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        Returns
        -------

        peaks : list
            A list containing the profile scores for each column of the given
            alignment.

        Examples
        --------


        """
        return [
            calign.score_profile(
                [k[i] for k in self._alm_matrix],
                [k[i] for k in self._alm_matrix],
                self.scorer,
                gap_weight=gap_weight)
            for i in range(len(self._alm_matrix[0]))
        ]

    def get_local_peaks(self, threshold=2, gap_weight=0.0):
        """
        Return all peaks in a given alignment.

        Parameters
        ----------
        threshold : { int, float } (default=2)
            The threshold to determine whether a given column is a peak or not.
        gap_weight : float (default=0.0)
            The weight for gaps.

        """
        peaks = self.get_peaks(gap_weight=gap_weight)
        self.local = [i for i in range(len(peaks)) if peaks[i] > threshold]

    def get_pairwise_alignments(self, **keywords):
        """
        Function creates a dictionary of all pairwise alignments  scores.

        Parameters
        ----------
        new_calc : bool (default=True)
            Specify, whether the analysis should be repeated from the
            beginning, or whether already conducted analyses should be carried
            out.

        model : string (default="sca")
            A string indicating the name of the :py:class:`Model \
            <lingpy.data.model>` object that shall be used for the analysis.
            Currently, three models are supported:

            * "dolgo" -- a sound-class model based on :evobib:`Dolgopolsky1986`,

            * "sca" -- an extension of the "dolgo" sound-class model based on
              :evobib:`List2012b`, and

            * "asjp" -- an independent sound-class model which is based on the
              sound-class model of :evobib:`Brown2008` and the empirical data
              of :evobib:`Brown2011` (see the description in
              :evobib:`List2012`.

        mode : string (default="global")
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

        gop : int (default=-3)
            The gap opening penalty (GOP) used in the analysis.

        gep_scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1982`.

        factor : float (default=1)
            The factor by which the initial and the descending position shall
            be modified.

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        restricted_chars : string (default="T")
            Define which characters of the prosodic string of a sequence
            reflect its secondary structure (cf. :evobib:`List2012b`) and
            should therefore be aligned specifically. This defaults to "T",
            since this is the character that represents tones in the prosodic
            strings of sequences.

        """
        setdefaults(
            keywords,
            new_calc=True,
            model=rcParams['sca'],
            mode='global',
            gop=-3,
            scale=0.5,
            factor=1,
            restricted_chars='T_',
            classes=True,
            sonar=True,
            scorer={})

        if keywords['new_calc']:
            # define the class model
            self._set_model(
                keywords['model'],
                keywords['classes'],
                keywords['sonar'],
                keywords['scorer'])

            # reset the scorer to "classes"
            self._set_scorer('classes')

            # retrieve the alignments
            self._get_pairwise_alignments(
                keywords['mode'],
                keywords['gop'],
                keywords['scale'],
                keywords['factor'],
                keywords['restricted_chars'])

            self.alignments = {}

            for i, j in combinations_with_replacement(range(self.height), 2):
                # get the score of the alignment
                score = self.matrix[i][j]

                # retrieve the numeric tokens
                tokA = self._alignments[i][j][0]
                tokB = self._alignments[i][j][1]

                # append values to dictionary
                for k in self.int2ext[i]:
                    for l in self.int2ext[j]:
                        almA, almB = [], []
                        for idx, tok, alm in [(k, tokA, almA), (l, tokB, almB)]:
                            for m in tok:
                                try:
                                    alm.append(self._get(
                                        dotjoin(idx + 1, m.split('.')[1]), 'tokens'))
                                except:
                                    alm.append('-')

                        self.alignments[k, l] = [almA, almB, score]
        else:
            # if new_calc is not chosen, the PID of an alignment will be
            # returned, beware only to calculate the pid for unique sequences
            # in order to save time and memory
            self.alignments = {}
            for seqA, seqB in combinations_with_replacement(self.uniseqs.keys(), 2):
                # get the score of the alignment
                almA = self.alm_matrix[self.uniseqs[seqA][0]]
                almB = self.alm_matrix[self.uniseqs[seqB][0]]

                score = pid(almA, almB, 2)

                for k, l in product(self.uniseqs[seqA], self.uniseqs[seqB]):
                    if k != l:
                        self.alignments[k, l] = [almA, almB, score]

    def get_pid(self, mode=1):
        """
        Return the Percentage Identity (PID) score of the calculated MSA.

        Parameters
        ----------
        mode : { 1, 2, 3, 4, 5 } (default=1)
            Indicate which of the four possible PID scores described in
            :evobib:`Raghava2006` should be calculated, the fifth possibility is added
            for linguistic purposes:

            1. identical positions / (aligned positions + internal gap positions),

            2. identical positions / aligned positions,

            3. identical positions / shortest sequence, or

            4. identical positions / shortest sequence (including internal gap
               pos.)

            5. identical positions / (aligned positions + 2 * number of gaps)


        Returns
        -------
        score : float
            The PID score of the given alignment as a floating point number between
            0 and 1.

        See also
        --------
        lingpy.sequence.sound_classes.pid

        """

        # create a dictionary of unique sequences
        weights = {}
        for key in self.uniseqs.keys():
            vals = self.uniseqs[key]
            l = len(vals)
            for val in vals:
                weights[val] = (key, 1.0 / l)

        score = 0.0
        for (i, seqA), (j, seqB) in combinations(enumerate(self.alm_matrix), 2):
            kA, wA = weights[i]
            kB, wB = weights[j]
            if kA != kB:
                score += pid(seqA, seqB, mode) * wA * wB

        l = len(self.uniseqs)
        count = (l ** 2 - l) / 2
        return 1.0 if count == 0 else score / count

    def _similar_gap_sites(self):
        """
        Create a dictionary which lists all strings with a similar gap
        structure and their index. The gap structure is marked in the key of
        the dictionary, where '1' refers to gapped sites and '0' refers to
        ungapped sites. The values of the dictionary are a list of integers
        referring to the position of the sequences having this structure in the
        MSA.
        """

        self.gap_dict = defaultdict(list)
        """
        A dictionary storing the different gap profiles of an MSA as keys and
        the indices of the corresponding sequences as values.
        """
        for i in range(len(self._classes)):
            gaps = ''.join('1' if seg == 'X' else '0' for seg in self._alm_matrix[i])
            self.gap_dict[gaps].append(i)

    def _mk_gap_array(self):
        """
        Return an array of the gap-profiles of an alignment.
        @return: An array, representing the gap profiles as integers (0
            indicates characters and 1 indicates gaps).
        @rtype: C{scipy.array}
        """
        self._similar_gap_sites()
        gap_array = [0 for _ in self._alm_matrix]

        for key in self.gap_dict:
            for i in self.gap_dict[key]:
                gap_array[i] = [int(x) for x in key]

        return gap_array

    def _swap_condition(self):
        """
        The condition for swaps to possibly occur in the alignment. These are
        the complementary sites in the alignment, which are extracted from the
        gap array.
        """
        swaps = []
        gap_array = self._mk_gap_array()
        for i in range(len(gap_array[0])):
            try:
                if 0 not in [a + b for a, b in zip(
                    [line[i] for line in gap_array], [line[i + 2] for line in gap_array]
                )]:
                    swaps.append((i))
            except:
                pass

        return swaps

    def _swap_check(self, ind, gap_weight=1.0, swap_penalty=-5, db=False):
        """
        Carry out a check for swapped regions.
        """
        # [i] We define two versions of the possibly swapped region, a first
        # ... one, where the original alignment is unchanged, and a second one,
        # ... where the alignment is shifted, i.e. the gaps are switched.

        matA = [[c for c in line] for line in self._alm_matrix]
        matB = [[c for c in line] for line in self._alm_matrix]

        for i in range(len(self._classes)):
            # [i] shift the gap of the first and third matrix
            if matA[i][ind] != 'X' and matA[i][ind + 2] == 'X':
                matA[i][ind + 2] = matA[i][ind + 1]
                matA[i][ind + 1] = matA[i][ind]
                matA[i][ind] = 'X'

            elif matA[i][ind] == 'X':
                matA[i][ind] = matA[i][ind + 1]
                matA[i][ind + 1] = matA[i][ind + 2]
                matA[i][ind + 2] = 'X'

        # determine in which direction to turn by counting the number of chars
        # in all cols
        t1 = len([i for i in [line[ind] for line in matA] if i != 'X'])
        t2 = len([i for i in [line[ind + 2] for line in matA] if i != 'X'])
        turnLeftA = t1 > t2

        t1 = len([i for i in [line[ind] for line in matB] if i != 'X'])
        t2 = len([i for i in [line[ind + 2] for line in matB] if i != 'X'])
        turnLeftB = t1 > t2

        for i in range(len(self._classes)):
            # [i] unswap the possibly swapped columns by shifting values
            # ... unequal to a gap and leaving a special symbol (+) which will
            # ... cope for the penalty for a swap.
            if turnLeftA:
                if matA[i][ind] != 'X':
                    pass
                elif matA[i][ind + 2] != 'X':
                    matA[i][ind] = matA[i][ind + 2]
                    matA[i][ind + 2] = '+'
            else:
                if matA[i][ind] != 'X':
                    matA[i][ind + 2] = matA[i][ind]
                    matA[i][ind] = '+'
                else:
                    pass

            # [i] apply the same procedure to the unshifted matrix
            if turnLeftB:
                if matB[i][ind] != 'X':
                    pass
                elif matB[i][ind + 2] != 'X':
                    matB[i][ind] = matB[i][ind + 2]
                    matB[i][ind + 2] = '+'
            else:
                if matB[i][ind] != 'X':
                    matB[i][ind + 2] = matB[i][ind]
                    matB[i][ind] = '+'
                else:
                    pass

        # [i] calculate normal and new sum-of-pairs scores, convert to integers
        # ... in order to guarantee the accuracy of the comparison of
        # ... sop-scores
        msa = int(self.sum_of_pairs(gap_weight=gap_weight) * 1000000)
        msaA = self._swap_sum_of_pairs(
            matA, gap_weight=gap_weight, swap_penalty=swap_penalty)
        msaB = self._swap_sum_of_pairs(
            matB, gap_weight=gap_weight, swap_penalty=swap_penalty)
        msaAB = int((msaA + msaB) * 0.5 * 1000000)
        self.log.debug('', extra=dict(
            lines=[[self._get(x, '_classes') for x in line] for line in matA] + \
                    [msaA] + \
                    [[self._get(x, '_classes') for x in line] for line in matA] +\
                    [msaB] +
                    [msaAB, msa]
            ))

        # return True if the newly calculated sop-score is greater than the previous one
        return msaAB > msa

    def swap_check(self, swap_penalty=-3, score_mode='classes'):
        """
        Check for possibly swapped sites in the alignment.

        Parameters
        ----------
        swap_penalty : { int, float } (default=-3)
            Specify the penalty for swaps in the alignment.

        score_mode : { "classes", "library" } (default="classes")
            Define the score-mode of the calculation which is either based on
            sound classes proper, or on the specific scores derived from the
            library approach.

        Returns
        -------
        result : bool
            Returns ``True``, if a swap was identified, and ``False``
            otherwise. The information regarding the position of the swap is
            stored in the attribute ``swap_index``.

        Notes
        -----
        The method for swap detection is described in detail in :evobib:`List2012b`.

        Examples
        --------
        Define a set of strings whose alignment contans a swap.

        >>> from lingpy import *
        >>> mult = Multiple(["woldemort", "waldemar", "wladimir"])

        Align the data, using the progressive approach.

        >>> mult.prog_align()

        Check for swaps.

        >>> mult.swap_check()
        True

        Print the alignment

        >>> print(mult)
        w   o   l   -   d   e   m   o   r   t
        w   a   l   -   d   e   m   a   r   -
        v   -   l   a   d   i   m   i   r   -

        """
        self._set_scorer(score_mode)
        swaps = []

        for i in self._swap_condition():
            if self._swap_check(i, swap_penalty=swap_penalty):
                swaps.append(i)
            else:
                pass

        if swaps:
            # check for incompatible swaps. this is a temporary solution
            # since it might be better to check the BEST swaps and not only
            # to discard them in a linear order
            swp = []
            while swaps:
                i = swaps.pop(0)
                if i - 1 not in swp and i - 2 not in swp and i - 3 not in swp:
                    swp.append(i)

            self.swap_index = [(i_, i_ + 1, i_ + 2) for i_ in swp]
            self.swaps = self.swap_index
            return True
        return False


def mult_align(seqs, gop=-1, scale=0.5, tree_calc='upgma', scoredict=False,
        pprint=False):
    """
    A short-cut method for multiple alignment analyses.

    Parameters
    ----------
    seqs : list
        The input sequences.
    gop = int (default=-1)
        The gap opening penalty.
    scale : float (default=0.5)
        The scaling factor by which penalties for gap extensions are decreased.
    tree_calc : { "upgma" "neighbor" } (default="upgma")
        The algorithm which is used for the calculation of the guide tree.
    pprint : bool (default=False)
        Indicate whether results shall be printed onto screen.

    Returns
    -------
    alignments : list
        A two-dimensional list in which alignments are represented as a list of
        tokens.

    Examples
    --------
    >>> m = mult_align(["woldemort", "waldemar", "vladimir"], pprint=True)
    w   o   l   -   d   e   m   o   r   t
    w   a   l   -   d   e   m   a   r   -
    -   v   l   a   d   i   m   i   r   -


    """
    m = Multiple(seqs)
    m.prog_align(
        classes=False,
        sonar=False,
        gop=gop,
        tree_calc=tree_calc,
        scale=scale,
        scoredict=scoredict or {})
    
    as_string(m, pprint=pprint)
    return m.alm_matrix
