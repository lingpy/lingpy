# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
import random
from itertools import combinations_with_replacement, product, combinations

from collections import Counter, defaultdict
from math import factorial
from copy import copy

# thirdparty
from six import text_type
from six import string_types
import numpy as np

# lingpy-modules
from ..settings import rcParams
from ..sequence.sound_classes import (
    ipa2tokens, tokens2class, prosodic_string, prosodic_weights, class2tokens)
from ..sequence.generate import MCPhon
from ..basic import Wordlist
from ..align.pairwise import turchin, edit_dist
from ..convert.strings import scorer2str
from ..algorithm import clustering
from ..algorithm import calign
from ..algorithm import talign
from ..algorithm import misc
from .. import util
from .. import log


def _check_tokens(key_and_tokens):
    """Generator for error reports on token strings.

    :param key_and_tokens: iterator over (key, token_string) pairs.
    """
    for key, line in key_and_tokens:
        if "" in line:
            yield (key, "empty token", line)
        else:
            try:
                sonars = tokens2class(line, rcParams['art'])
                if not sonars or sonars == ['0']:
                    yield (key, "empty sound-class string", line)
                elif '0' in sonars:
                    yield (key, "bad character in tokens at «{0}»".format(
                        line[sonars.index('0')]), line)
            except ValueError:
                yield (key, "sound-class conversion failed", line)


def charstring(id_, char='X', cls='-'):
    return '{0}.{1}.{2}'.format(id_, char, cls)


def char_from_charstring(cstring):
    comps = cstring.split('.')
    if len(comps) == 3:
        # a full charstring
        return comps[1][0]
    if len(comps) == 2:
        # a reduced charstring
        return comps[0][0]
    raise ValueError(cstring)


def get_score_dict(chars, model):
    matrix = [[0.0 for i in chars] for j in chars]
    for (i, charA), (j, charB) in combinations_with_replacement(
            enumerate(chars), r=2):
        matrix[i][j] = model(char_from_charstring(charA), char_from_charstring(charB))
        if i < j:
            matrix[j][i] = matrix[i][j]
    return misc.ScoreDict(chars, matrix)


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

    Notes
    -----
    Instantiating this class does not require a lot of parameters. However,
    the user may modify its behaviour by providing additional attributes in the
    input file.

    """
    def __init__(self, filename, **keywords):
        kw = {
            "model": rcParams['sca'],
            "merge_vowels": rcParams['merge_vowels'],
            'transform': rcParams['lexstat_transform'],
            "check": False,
            "apply_checks": False,
            "defaults": False,
            "no_bscorer": False,
            "errors": "errors.log",
        }
        kw.update(keywords)

        if isinstance(kw['model'], string_types):
            self.model = rcParams[kw['model']]
        else:
            self.model = kw['model']

        # set the lexstat stamp
        self._stamp = "# Created using the LexStat class of LingPy-2.0\n"

        # initialize the wordlist
        Wordlist.__init__(self, filename)
        assert "tokens" in self.header or "ipa" in self.header

        # check for basic input data
        # tokens
        if "tokens" not in self.header:
            self.add_entries(
                "tokens", "ipa", lambda x: ipa2tokens(x, merge_vowels=kw['merge_vowels']))

        # add a debug procedure for tokens
        if kw["check"]:
            errors = list(_check_tokens((key, self[key, "tokens"]) for key in self))
            if errors:
                lines = ["ID\tTokens\tError-Type"]
                for key, msg, line in errors:
                    lines.append("{0}\t<{1}>\t{2}".format(key, msg, ' '.join(line)))
                util.write_text_file(kw['errors'], lines)

                if kw["apply_checks"] or util.confirm(
                        "There were errors in the input data - exclude them?"):
                    self.output(
                        'qlc',
                        filename=self.filename + '_cleaned',
                        subset=True,
                        rows={"ID": "not in " + str([i[0] for i in errors])})
                    # load the data in a new LexStat instance and copy the __dict__
                    lexstat = LexStat(self.filename + '_cleaned.qlc', **kw)
                    lexstat._meta['errors'] = [i[0] for i in errors]
                    self.__dict__ = copy(lexstat.__dict__)
                return
            else:
                log.info("No obvious errors found in the data.")

        # sonority profiles
        if "sonars" not in self.header:
            self.add_entries(
                "sonars",
                "tokens",
                lambda x: [int(i) for i in tokens2class(
                    x, rcParams['art'], stress=rcParams['stress'])])
        if "prostrings" not in self.header:
            self.add_entries("prostrings", "sonars", lambda x: prosodic_string(x))
        # get sound class strings
        if "classes" not in self.header:
            self.add_entries(
                "classes", "tokens", lambda x: ''.join(tokens2class(x, kw["model"])))
        # create IDs for the languages
        if "langid" not in self.header:
            transform = dict(zip(self.taxa, [str(i + 1) for i in range(self.width)]))
            self.add_entries("langid", "taxa", lambda x: transform[x])
        # get the numbers for all strings
        if "numbers" not in self.header:
            # change the discriminative potential of the sound-class string
            # tuples, note that this is still wip, we have to tweak around with
            # this in order to find an optimum for the calculation
            self._transform = kw['transform']
            self.add_entries(
                "numbers",
                "langid,classes,prostrings",
                lambda x, y: [charstring(x[y[0]], a, self._transform[b])
                              for a, b in zip(x[y[1]], x[y[2]])])
        # check for weights
        if "weights" not in self.header:
            self.add_entries("weights", "prostrings", lambda x: prosodic_weights(x))

        # check for duplicates
        # first, check for item 'words' in data, if this is not given, create it
        if 'ipa' not in self.header:
            self.add_entries('ipa', 'tokens', lambda x: ''.join(x))

        if "duplicates" not in self.header:
            duplicates = {}
            for taxon in self.taxa:
                words = set()
                for idx in self.get_list(col=taxon, flat=True):
                    word = self[idx, 'ipa']
                    duplicates[idx] = 1 if word in words else 0
                    words.add(word)
            self.add_entries("duplicates", duplicates, lambda x: x)

        # create an index
        if not hasattr(self, 'freqs'):
            self.chars = set()
            self.freqs = {}

            for taxon in self.taxa:
                self.freqs[taxon] = Counter()
                for word in self.get_list(col=taxon, entry='numbers', flat=True):
                    self.freqs[taxon].update(word)
                self.chars = self.chars.union(self.freqs[taxon].keys())

            self.rchars = sorted(set(char.split('.', 1)[1] for char in self.chars))
            self.chars = sorted(self.chars) \
                + [charstring(i + 1) for i in range(self.width)]

        if not hasattr(self, "scorer"):
            self._meta['scorer'] = {}

        # create a scoring dictionary
        if not hasattr(self, "bscorer") and not kw['no_bscorer']:
            self._meta['scorer']['bscorer'] = self.bscorer = get_score_dict(
                self.chars, self.model)
        elif not kw['no_bscorer']:
            self.bscorer = self._meta['scorer']['bscorer']

        if not hasattr(self, "rscorer"):
            self._meta['scorer']['rscorer'] = self.rscorer = get_score_dict(
                self.rchars, self.model)

        if 'scorer' in self._meta:
            if 'cscorer' in self._meta['scorer']:
                self.cscorer = self._meta['scorer']['cscorer']

        # make the language pairs
        if not hasattr(self, "pairs"):
            self.pairs = {}
            for (i, taxonA), (j, taxonB) in combinations_with_replacement(
                    enumerate(self.taxa), r=2):
                self.pairs[taxonA, taxonB] = []
                dictA = self.get_dict(col=taxonA)
                dictB = self.get_dict(col=taxonB)
                if i < j:
                    for c in sorted(set(dictA).intersection(dictB)):
                        for idxA, idxB in product(dictA[c], dictB[c]):
                            dA = self[idxA, "duplicates"]
                            dB = self[idxB, "duplicates"]
                            if dA != 1 and dB != 1:
                                self.pairs[taxonA, taxonB] += [(idxA, idxB)]
                elif i == j:
                    for c in sorted(dictA):
                        for idx in dictA[c]:
                            dAB = self[idx, "duplicates"]
                            if dAB != 1:
                                self.pairs[taxonA, taxonA] += [(idx, idx)]

    def __repr__(self):
        return "<lexstat-model {0}>".format(self.filename)

    def __getitem__(self, idx):
        """
        Method allows quick access to the data by passing the integer key.

        In contrast to the basic wordlist, the LexStat wordlist further allows
        to access item pairs by passing a tuple.
        """
        if idx in self._cache:
            return self._cache[idx]

        if idx in self._data:
            self._cache[idx] = self._data[idx]
            return self._cache[idx]

        try:
            return (
                self._data[idx[0][0]][self._header[self._alias[idx[1]]]],
                self._data[idx[0][1]][self._header[self._alias[idx[1]]]])
        except:
            try:
                # return data entry with specified key word
                self._cache[idx] = self._data[idx[0]][self._header[self._alias[idx[1]]]]
                return self._cache[idx]
            except KeyError:
                pass

    def get_subset(self, sublist, ref='concept'):
        """
        Function creates a specific subset of all word pairs.

        Parameters
        ----------
        sublist : list
            A list which contains those items which should be considered for
            the subset creation, for example, a list of concepts.
        ref : string (default="concept")
            The reference point to compare the given sublist.

        Notes
        -----
        This function can be used to consider only a smaller part of word pairs
        when creating a scorer. Normally, all words are compared, but defining
        a subset allows to compare only those belonging to a specific concept
        list (Swadesh list).
        """
        self.subsets = {}
        for tA, tB in combinations_with_replacement(self.taxa, r=2):
            self.subsets[tA, tB] = [
                pair for pair in self.pairs[tA, tB] if self[pair, ref][0] in sublist]

    def _get_corrdist(self, **keywords):
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
            preprocessing_threshold=rcParams['lexstat_preprocessing_threshold'],
            ref='scaid',
            restricted_chars=rcParams['restricted_chars'],
            threshold=rcParams['lexstat_threshold'],
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

        tasks = factorial(len(self.taxa) + 1) / 2 / factorial(len(self.taxa) - 1)
        with util.ProgressBar('CORRESPONDENCE CALCULATION', tasks) as progress:
            for (i, tA), (j, tB) in combinations_with_replacement(
                    enumerate(self.taxa), r=2):
                progress.update()
                log.info("Calculating alignments for pair {0} / {1}.".format(tA, tB))

                corrdist[tA, tB] = defaultdict(float)
                for mode, gop, scale in kw['modes']:
                    # XXX this is where we should add the new function for
                    # subsets of swadesh lists XXX
                    # this can be easily done by first checking for a
                    # sublist parameter and then getting all the numbers in
                    # a temporary variable "pairs" for all cases where this
                    # subset is defined, all that needs to be done is to
                    # provide an extra function that creates a
                    # subset-variable or hash in which for all language
                    # pairs the subset is defined.
                    pairs = self.pairs[tA, tB]
                    if kw['subset']:
                        pairs = [pair for pair in pairs if pair in self.subsets[tA, tB]]

                    if kw['preprocessing']:
                        pairs = [pair for pair in pairs
                                 if self[pair, kw['ref']][0] == self[pair, kw['ref']][1]]
                        threshold = 10.0
                    else:
                        threshold = kw['preprocessing_threshold']

                    corrs, self._included[tA, tB] = calign.corrdist(
                        threshold,
                        [self[pair, "numbers"] for pair in pairs],
                        [self[pair, "weights"] for pair in pairs],
                        [self[pair, "prostrings"] for pair in pairs],
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
                            a = charstring(i + 1)
                        elif b == '-':
                            b = charstring(j + 1)
                        corrdist[tA, tB][a, b] += d / len(kw['modes'])

        return corrdist

    def _get_randist(self, **keywords):
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
        tasks = factorial(len(self.taxa) + 1) / 2 / factorial(len(self.taxa) - 1)

        if method == 'markov':
            seqs, pros, weights = {}, {}, {}

            # get a random distribution for all pairs
            sample = random.sample(
                [(i, j) for i in range(kw['rands']) for j in range(kw['rands'])],
                kw['runs'])

            with util.ProgressBar('SEQUENCE GENERATION', len(self.taxa)) as progress:
                for i, taxon in enumerate(self.taxa):
                    progress.update()
                    log.info("Analyzing taxon {0}.".format(taxon))

                    tokens = self.get_list(col=taxon, entry="tokens", flat=True)
                    prostrings = self.get_list(col=taxon, entry="prostrings", flat=True)
                    m = MCPhon(tokens, True, prostrings)
                    words = []
                    j, k = 0, 0
                    while j < kw['rands']:
                        s = m.get_string(new=False)
                        if s in words:
                            k += 1
                            if k > kw['limit']:
                                break
                        else:
                            j += 1
                            words += [s]

                    seqs[taxon], pros[taxon], weights[taxon] = [], [], []
                    for w in words:
                        cls = tokens2class(w.split(' '), self.model)
                        pros[taxon].append(prosodic_string(w.split(' ')))
                        weights[taxon].append(prosodic_weights(pros[taxon][-1]))
                        seqs[taxon].append(
                            ['{0}.{1}'.format(c, p) for c, p in
                             zip(cls, [self._transform[pr] for pr in pros[taxon][-1]])])

            with util.ProgressBar('RANDOM CORRESPONDENCE CALCULATION', tasks) as progress:
                for (i, tA), (j, tB) in combinations_with_replacement(
                        enumerate(self.taxa), r=2):
                    progress.update()
                    log.info(
                        "Calculating random alignments for pair {0} / {1}.".format(tA, tB)
                    )
                    corrdist[tA, tB] = defaultdict(float)
                    for mode, gop, scale in kw['modes']:
                        corrs, included = calign.corrdist(
                            10.0,
                            [(seqs[tA][x], seqs[tB][y]) for x, y in sample],
                            [(weights[tA][x], weights[tB][y]) for x, y in sample],
                            [(pros[tA][x], pros[tB][y]) for x, y in sample],
                            gop,
                            scale,
                            kw['factor'],
                            self.rscorer,
                            mode,
                            kw['restricted_chars'])

                        # change representation of gaps
                        for a, b in list(corrs.keys()):
                            # get the correspondence count
                            d = corrs[a, b] * self._included[tA, tB] / included
                            # XXX check XXX * len(self.pairs[tA,tB]) / runs

                            # check for gaps
                            if a == '-':
                                a = 'X.-'
                            elif b == '-':
                                b = 'X.-'

                            a = str(i + 1) + '.' + a
                            b = str(j + 1) + '.' + b
                            corrdist[tA, tB][a, b] += d / len(kw['modes'])
        # use shuffle approach otherwise
        else:
            tasks = factorial(len(self.taxa) + 1) / 2 / factorial(len(self.taxa) - 1)
            with util.ProgressBar('RANDOM CORRESPONDENCE CALCULATION', tasks) as progress:
                for (i, tA), (j, tB) in combinations_with_replacement(
                        enumerate(self.taxa), r=2):
                    progress.update()
                    log.info(
                        "Calculating random alignments for pair {0} / {1}.".format(tA, tB)
                    )
                    corrdist[tA, tB] = defaultdict(float)

                    # get the number pairs etc.
                    numbers = [self[pair, 'numbers'] for pair in self.pairs[tA, tB]]
                    gops = [self[pair, 'weights'] for pair in self.pairs[tA, tB]]
                    prostrings = [self[pair, 'prostrings'] for pair in self.pairs[tA, tB]]

                    sample = [(x, y)
                              for x in range(len(numbers)) for y in range(len(numbers))]
                    if len(sample) > kw['runs']:
                        sample = random.sample(sample, kw['runs'])

                    for mode, gop, scale in kw['modes']:
                        corrs, included = calign.corrdist(
                            10.0,
                            [(numbers[s[0]][0], numbers[s[1]][1]) for s in sample],
                            [(gops[s[0]][0], gops[s[1]][1]) for s in sample],
                            [(prostrings[s[0]][0], prostrings[s[1]][1]) for s in sample],
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
                                a = charstring(i + 1)
                            elif b == '-':
                                b = charstring(j + 1)

                            corrdist[tA, tB][a, b] += d / len(kw['modes'])
        return corrdist

    def get_scorer(self, **keywords):
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
        kw = dict(
            method=rcParams['lexstat_scoring_method'],
            ratio=rcParams['lexstat_ratio'],
            vscale=rcParams['lexstat_vscale'],
            runs=rcParams['lexstat_runs'],
            threshold=rcParams['lexstat_threshold'],
            modes=rcParams['lexstat_modes'],
            factor=rcParams['align_factor'],
            restricted_chars=rcParams['restricted_chars'],
            force=False,
            preprocessing=True,
            rands=rcParams['lexstat_rands'],
            limit=rcParams['lexstat_limit'],
            cluster_method=rcParams['lexstat_cluster_method'],
            gop=rcParams['align_gop'],
            preprocessing_threshold=rcParams['lexstat_preprocessing_threshold'],
            preprocessing_method=rcParams['lexstat_preprocessing_method'],
            subset=False,
            defaults=False,
        )
        kw.update(keywords)
        if kw['defaults']:
            return kw

        # get parameters and store them in string
        params = dict(
            ratio=kw['ratio'],
            vscale=kw['vscale'],
            runs=kw['runs'],
            threshold=kw['preprocessing_threshold'],
            modestring=':'.join(
                '{0}-{1}-{2:.2f}'.format(a, abs(b), c) for a, b, c in kw['modes']),
            factor=kw['factor'],
            restricted_chars=kw['restricted_chars'],
            method=kw['method'],
            preprocessing='{0}:{1}:{2}'.format(
                kw['preprocessing'], kw['cluster_method'], kw['gop']))

        parstring = '_'.join(
            [
                '{ratio[0]}:{ratio[1]}'
                '{vscale:.2f}',
                '{runs}',
                '{threshold:.2f}',
                '{modestring}',
                '{factor:.2f}',
                '{restricted_chars}',
                '{method}',
                '{preprocessing}'
            ]).format(**params)

        # check for existing attributes
        if hasattr(self, 'cscorer') and not kw['force']:
            log.warn(
                "An identical scoring function has already been calculated, force "
                "recalculation by setting 'force' to 'True'.")
            return

        # check for attribute
        if hasattr(self, 'params') and not kw['force']:
            if 'cscorer' in self.params:
                if self.params['cscorer'] == params:
                    log.warn(
                        "An identical scoring function has already been calculated, force"
                        " recalculation by setting 'force' to 'True'.")
                    return
            else:
                log.warn(
                    "A different scoring function has already been calculated, "
                    "overwriting previous settings.")

        # store parameters
        self.params = {'cscorer': params}
        self._meta['params'] = self.params
        self._stamp += "# Parameters: " + parstring + '\n'

        # get the correspondence distribution
        self._corrdist = self._get_corrdist(**kw)
        # get the random distribution
        self._randist = self._get_randist(**kw)

        # get the average gop
        gop = sum([m[1] for m in kw['modes']]) / len(kw['modes'])

        # create the new scoring matrix
        matrix = [[c for c in line] for line in self.bscorer.matrix]
        char_dict = self.bscorer.chars2int

        for (i, tA), (j, tB) in combinations_with_replacement(enumerate(self.taxa), r=2):
            for charA, charB in product(
                list(self.freqs[tA]) + [charstring(i + 1)],
                list(self.freqs[tB]) + [charstring(j + 1)]
            ):
                exp = self._randist.get((tA, tB), {}).get((charA, charB), False)
                att = self._corrdist.get((tA, tB), {}).get((charA, charB), False)

                # in the following we follow the former lexstat protocol
                if att <= 1 and i != j:
                    att = False

                if att and exp:
                    score = np.log2((att ** 2) / (exp ** 2))
                elif att and not exp:
                    score = np.log2((att ** 2) / 0.00001)
                elif exp and not att:
                    score = -5  # XXX gop ???
                else:  # elif not exp and not att:
                    score = -90  # ???

                # combine the scores
                if '-' not in charA + charB:
                    sim = self.bscorer[charA, charB]
                else:
                    sim = gop

                # get the real score
                rscore = (kw['ratio'][0] * score + kw['ratio'][1] * sim) \
                    / sum(kw['ratio'])

                try:
                    idxA = char_dict[charA]
                    idxB = char_dict[charB]

                    # use the vowel scale
                    if charA[4] in 'XYZT_' and charB[4] in 'XYZT_':
                        matrix[idxA][idxB] = matrix[idxB][idxA] = kw['vscale'] * rscore
                    else:
                        matrix[idxA][idxB] = matrix[idxB][idxA] = rscore
                except:
                    pass

        self.cscorer = misc.ScoreDict(self.chars, matrix)
        self._meta['scorer']['cscorer'] = self.cscorer

    def align_pairs(self, idxA, idxB, concept=None, **keywords):
        """
        Align all or some words of a given pair of languages.

        Parameters
        ----------
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
        kw = dict(
            method='lexstat',
            mode="global",
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            pprint=True,
            return_distance=False,
            gop=-2,
            distance=True,
            defaults=False,
            return_raw=False
        )
        kw.update(keywords)
        if kw['defaults']:
            return kw

        if isinstance(idxA, (text_type, tuple)):
            if isinstance(idxA, tuple):
                for indexA, indexB in product(
                    self.get_dict(col=idxA[0])[idxA[1]],
                    self.get_dict(col=idxB[0])[idxB[1]],
                ):
                    self.align_pairs(indexA, indexB, **kw)
            else:
                if not concept:
                    for c in self.concepts:
                        print("Concept: {0}".format(c))
                        self.align_pairs(idxA, idxB, c, **kw)
                        print('')
                else:
                    self.align_pairs((idxA, concept), (idxB, concept), concept=None, **kw)
            return

        # assign the distance value
        distance = 1 if kw['distance'] else 0

        # get the language ids
        lA = self[idxA, 'langid']
        lB = self[idxB, 'langid']

        if kw['method'] == 'lexstat':
            scorer = self.cscorer
            gop = 1.0
            weightsA = [self.cscorer[charstring(lA), n] for n in self[idxA, 'numbers']]
            weightsB = [self.cscorer[charstring(lB), n] for n in self[idxB, 'numbers']]
        else:
            gop = kw['gop']
            weightsA = self[idxA, 'weights']
            weightsB = self[idxB, 'weights']
            scorer = self.bscorer

        almA, almB, d = calign.align_pair(
            self[idxA, 'numbers'],
            self[idxB, 'numbers'],
            weightsA,
            weightsB,
            self[idxA, 'prostrings'],
            self[idxB, 'prostrings'],
            gop,
            kw['scale'],
            kw['factor'],
            scorer,
            kw['mode'],
            kw['restricted_chars'],
            distance)

        # get a string of scores
        if kw['method'] == 'lexstat':
            fun = lambda x, y: x if x != '-' else charstring(y)
            scoreA = [fun(a, lA) for a in almA]
            scoreB = [fun(b, lB) for b in almB]
        else:
            scoreA = almA
            scoreB = almB

        scores = ['{0:.2f}'.format(scorer[a, b]) for a, b in zip(scoreA, scoreB)]

        if kw['return_raw']:
            return almA, almB, d

        almA = class2tokens(self[idxA, 'tokens'], almA)
        almB = class2tokens(self[idxB, 'tokens'], almB)
        if kw['pprint']:
            print('\t'.join(almA))
            print('\t'.join(almB))
            print('\t'.join(scores))
            if kw['distance']:
                print('Distance: {0:.2f}'.format(d))
            else:
                print('Similarity: {0:.2f}'.format(d))

        if kw['return_distance']:
            return d
        return almA, almB, d

    def _get_matrices(
            self,
            concept=False,
            method='sca',
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            mode='overlap',
            gop=-2,
            restriction='',
            **keywords):
        """
        Calculate alignment matrices.

        Notes
        -----
        This is an iterator object and it yields the indices of a given
        concept, the matrix, and the concept.
        """
        # currently, there are no defaults XXX
        kw = dict(
            defaults=False,
            external_scorer=False,  # external scoring function
        )
        kw.update(keywords)

        # check for method
        if method == 'lexstat':
            # check for scorer
            if not hasattr(self, 'cscorer'):
                log.warn("No correspondence-scorer has been specified.")
                return

            # define the function with help of lambda
            function = lambda idxA, idxy: calign.align_pair(
                self[idxA, 'numbers'],
                self[idxB, 'numbers'],
                [self.cscorer[charstring(self[idxB, 'langid']), n]
                 for n in self[idxA, 'numbers']],
                [self.cscorer[charstring(self[idxA, 'langid']), n]
                 for n in self[idxB, 'numbers']],
                self[idxA, 'prostrings'],
                self[idxB, 'prostrings'],
                1,
                scale,
                factor,
                self.cscorer,
                mode,
                restricted_chars,
                1)[2]
        elif method == 'sca':
            # define the function with help of lambda
            function = lambda idxA, idxB: calign.align_pair(
                [n.split('.', 1)[1] for n in self[idxA, 'numbers']],
                [n.split('.', 1)[1] for n in self[idxB, 'numbers']],
                self[idxA, 'weights'],
                self[idxB, 'weights'],
                self[idxA, 'prostrings'],
                self[idxB, 'prostrings'],
                gop,
                scale,
                factor,
                self.rscorer,
                mode,
                restricted_chars,
                1)[2]
        elif method == 'edit-dist':
            entry = kw.get('entry', 'tokens')
            function = lambda idxA, idxB: edit_dist(
                self[idxA, entry], self[idxB, entry], True, restriction)
        elif method == 'turchin':
            function = lambda idxA, idxB: turchin(
                self[idxA, 'tokens'], self[idxB, 'tokens'])
        elif method == 'custom':
            function = lambda idxA, idxB: talign.align_pair(
                self[idxA, 'utokens'],
                self[idxB, 'utokens'],
                gop,
                scale,
                keywords['external_scorer'],
                'overlap',
                True)[2]
        else:  # pragma: no cover
            raise ValueError(method)

        concepts = [concept] if concept else sorted(self.rows)
        for c in concepts:
            log.info("Analyzing words for concept <{0}>.".format(c))
            indices = self.get_list(row=c, flat=True)
            matrix = []
            for idxA, idxB in combinations(indices, r=2):
                try:
                    d = function(idxA, idxB)
                except ZeroDivisionError:
                    log.warn(
                        "Encountered Zero-Division for the comparison of "
                        "{0} and {1}".format(
                            ''.join(self[idxA, "tokens"]),
                            ''.join(self[idxB, "tokens"])))
                    d = 100

                matrix += [d]

            matrix = misc.squareform(matrix)

            if not concept:
                yield c, indices, matrix
            else:
                yield matrix

    def cluster(
            self,
            method='sca',
            cluster_method='upgma',
            threshold=0.3,
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            mode='overlap',
            gop=-2,
            restriction='',
            ref='',
            external_function=None,
            **keywords):
        """
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
            inflation=2,
            expansion=2,
            max_steps=1000,
            add_self_loops=True,
            guess_threshold=False,
            gt_trange=(0.4, 0.6, 0.02),
            mcl_logs=lambda x: -np.log2((1 - x) ** 2),
            gt_mode='average',
            matrix_type='distances',
            link_threshold=False,
            _return_matrix=False,  # help function for test purposes
            defaults=False,
            external_scorer=False,  # external scoring dictionary
        )
        kw.update(keywords)
        if kw['defaults']:
            return kw

        # check for parameters and add clustering, in order to make sure that
        # analyses are not repeated
        if not hasattr(self, 'params'):
            self.params = {}

        self.params['cluster'] = "{0}_{1}_{2:.2f}".format(
            method, cluster_method, threshold)
        self._stamp += '# Cluster: ' + self.params['cluster']

        if method not in ['lexstat', 'sca', 'turchin', 'edit-dist', 'custom']:
            raise ValueError("[!] The method you selected is not available.")

        # set up clustering algorithm, first the simple basics
        if external_function:
            fclust = external_function
        elif cluster_method in ['upgma', 'single', 'complete', 'ward']:
            fclust = lambda x, y: clustering.flat_cluster(
                cluster_method, y, x, revert=True)
        elif cluster_method == 'mcl':
            fclust = lambda x, y: clustering.mcl(
                y,
                x,
                list(range(len(x))),
                max_steps=kw['max_steps'],
                inflation=kw['inflation'],
                expansion=kw['expansion'],
                add_self_loops=kw['add_self_loops'],
                logs=kw['mcl_logs'],
                revert=True)
        elif cluster_method in ['lcl', 'link_clustering', 'lc']:
            fclust = lambda x, y: clustering.link_clustering(
                y,
                x,
                list(range(len(x))),
                revert=True,
                fuzzy=False,
                matrix_type=kw['matrix_type'],
                link_threshold=kw['link_threshold'])

        # make a dictionary that stores the clusters for later update
        clr = {}
        k = 0

        # create a matrix iterator
        matrices = self._get_matrices(
            method=method,
            scale=scale,
            factor=factor,
            restricted_chars=restricted_chars,
            mode=mode,
            gop=gop,
            restriction=restriction,
            **kw)

        if kw['guess_threshold']:
            thresholds = []
            # check for full consideration of basic t
            if kw['gt_mode'] == 'average':
                matrices = list(matrices)
                for c, i, m in matrices:
                    thresholds.append(clustering.best_threshold(m, kw['gt_trange']))

            # new method for threshold estimation based on calculating approximate
            # random distributions of similarities for each sequence
            elif kw['gt_mode'] == 'nulld':
                align = lambda x, y: self.align_pairs(
                    x,
                    y,
                    method=method,
                    restricted_chars=restricted_chars,
                    mode=mode,
                    scale=scale,
                    factor=factor,
                    return_distance=True,
                    pprint=False,
                    gop=gop)
                for l1, l2 in self.pairs:
                    if l1 != l2:
                        pairs = self.pairs[l1, l2]
                        for p1, p2 in pairs:
                            dx = [align(p1, pairs[random.randint(0, len(pairs) - 1)][1])
                                  for i in range(len(pairs) // 5)]
                            thresholds.extend(dx)
            if thresholds:
                threshold = sum(thresholds) / len(thresholds)

        with util.ProgressBar('SEQUENCE CLUSTERING', len(self.rows)) as progress:
            for concept, indices, matrix in matrices:
                progress.update()

                # check for keyword to guess the threshold
                if kw['guess_threshold'] and kw['gt_mode'] == 'item':
                    t = clustering.best_threshold(matrix, kw['gt_trange'])
                # FIXME: considering new function here JML
                #elif kw['guess_threshold'] and kw['gt_mode'] == 'nullditem':
                #    pass
                else:
                    t = threshold

                c = fclust(matrix, t)

                # specific clustering for fuzzy methods, currently not yet
                # supported
                #if cluster_method in ['fuzzy']:  # ['link_communities','lc','lcl']:
                #    clusters = [[d + k for d in c[i]] for i in range(len(matrix))]
                #    tests = []
                #    for clrx in clusters:
                #        for x in clrx:
                #            tests += [x]
                #    k = max(tests)
                #    for idxA, idxB in zip(indices, clusters):
                #        clr[idxA] = idxB
                #else:
                if 1:
                    # extract the clusters
                    clusters = [c[i] + k for i in range(len(matrix))]

                    # reassign the "k" value
                    k = max(clusters)

                    # add values to cluster dictionary
                    for idxA, idxB in zip(indices, clusters):
                        clr[idxA] = idxB

        if not ref:
            ref = method + 'id' if method in ['turchin', 'lexstat', 'sca', 'custom'] \
                else 'editid'
        self.add_entries(ref, clr, util.identity, override=kw.get('override', False))

        # assign thresholds to parameters
        self._current_threshold = threshold

    def _get_distances(
            self, method, mode, scale, factor, gop, sample, edit_dist_normalized):
        """
        :param sample: Callable returning an iterator of pairs sampled from the list of \
        pairs passed as sole argument.
        :param edit_dist_normalized: Whether edit_dist should be normalized.
        :return: generator of lists of distances for sampled pairs per taxa pair.
        """
        if method in ['sca', 'lexstat']:
            function = lambda x, y: self.align_pairs(
                x,
                y,
                method=method,
                distance=True,
                return_distance=True,
                pprint=False,
                mode=mode,
                scale=scale,
                factor=factor,
                gop=gop)
        else:
            function = lambda x, y: edit_dist(
                self[x, 'tokens'], self[y, 'tokens'], normalized=edit_dist_normalized)

        for taxA, taxB in combinations(self.taxa, r=2):
            distances = []
            for pA, pB in sample(self.pairs[taxA, taxB]):
                try:
                    d = function(pA, pB)
                except ZeroDivisionError:
                    self.log.error("Zero-Warning")
                    d = 1.0
                distances.append(d)
            yield distances

    def get_random_distances(
            self,
            method='lexstat',
            runs=100,
            mode='overlap',
            gop=-2,
            scale=0.5,
            factor=0.3,
            restricted_chars='T_'):
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
        def sample(pairs):
            _sample = random.sample(
                [(x, y) for x in range(len(pairs)) for y in range(len(pairs))],
                len(pairs))
            if len(_sample) > runs:
                _sample = random.sample(_sample, runs)
            return [(pairs[x][0], pairs[y][1]) for x, y in _sample]

        D = []
        for distances in self._get_distances(
                method, mode, scale, factor, gop, sample, False):
            D.extend(distances)
        return sorted(D)

    def get_distances(
            self,
            method='sca',
            mode='overlap',
            gop=-2,
            scale=0.5,
            factor=0.3,
            restricted_chars='T_',
            aggregate=True):
        """
        Method calculates different distance estimates for language pairs.

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
        aggregate : bool (default=True)
            Return aggregated distances in form of a distance matrix for all
            taxa in the data.

        Returns
        -------
        D : c{numpy.array}
            An array with all distances calculated for each sequence pair.
        """
        D = []
        for distances in self._get_distances(
                method, mode, scale, factor, gop, util.identity, True):
            if aggregate:
                D.append(sum(distances) / len(distances))
            else:
                D.extend(distances)
        return misc.squareform(D) if aggregate else sorted(D)

    def get_frequencies(self, ftype='sounds', ref='tokens', aggregated=False):
        """
        Computes the frequencies of a given wordlist.

        Parameters
        ----------
        ftype: str (default='sounds')
            The type of frequency which shall be calculated. Select between
            "sounds" (type-token frequencies of sounds), and "wordlength" (average
            word length per taxon or in aggregated form), or "diversity" for the diversity
            index (requires that you have carried out cognate judgments, and
            make sure to set the "ref" keyword to the column in which your
            cognates are).
        ref : str (default="tokens")
            The reference column, with the column for "tokens" as a default.
            Make sure to modify this keyword in case you want to check for the
            "diversity".
        aggregated : bool (default=False)
            Determine whether frequencies should be calculated in an aggregated
            way, for all languages, or on a language-per-language basis.

        Returns
        -------
        freqs : {dict, float}
            Depending on the selection of the datatype you chose, this returns
            either a dictionary containing the frequencies or a float
            indicating the ratio.
        """
        if ftype == 'sounds':
            _F = defaultdict(Counter)
            F = Counter()
            for k in self:
                tokens = self[k, ref]
                if tokens:
                    _F[self[k][self._colIdx]].update(tokens)
                    F.update(tokens)
            return F if aggregated else _F

        if ftype == 'wordlength':
            _W = defaultdict(lambda: [0, 0])
            for k in self:
                _W[self[k][self._colIdx]][0] += len(self[k, ref])
                _W[self[k][self._colIdx]][1] += 1
            _W = {a: b[0] / b[1] for a, b in _W.items()}
            return sum(_W.values()) / self.width if aggregated else _W

        if ftype == 'diversity':
            return (len(self.get_etymdict(ref)) - self.height) / (len(self) - self.height)

    def output(self, fileformat, **keywords):
        """
        Write data for lexstat to file.

        Parameters
        ----------
        fileformat : {'csv', 'tre','nwk','dst', 'taxa','starling', 'paps.nex', 'paps.csv'}
            The format that is written to file. This corresponds to the file
            extension, thus 'csv' creates a file in csv-format, 'dst' creates
            a file in Phylip-distance format, etc.
        filename : str
            Specify the name of the output file (defaults to a filename that
            indicates the creation date).
        subset : bool (default=False)
            If set to c{True}, return only a subset of the data. Which subset
            is specified in the keywords 'cols' and 'rows'.
        cols : list
            If *subset* is set to c{True}, specify the columns that shall be
            written to the csv-file.
        rows : dict
            If *subset* is set to c{True}, use a dictionary consisting of keys
            that specify a column and values that give a Python-statement in
            raw text, such as, e.g., "== 'hand'". The content of the specified
            column will then be checked against statement passed in the
            dictionary, and if it is evaluated to c{True}, the respective row
            will be written to file.
        cognates : str
            Name of the column that contains the cognate IDs if 'starling' is
            chosen as an output format.

        missing : { str, int } (default=0)
            If 'paps.nex' or 'paps.csv' is chosen as fileformat, this character
            will be inserted as an indicator of missing data.

        tree_calc : {'neighbor', 'upgma'}
            If no tree has been calculated and 'tre' or 'nwk' is chosen as
            output format, the method that is used to calculate the tree.

        threshold : float (default=0.6)
            The threshold that is used to carry out a flat cluster analysis if
            'groups' or 'cluster' is chosen as output format.
        """
        kw = dict(filename=self.filename, defaults=False)
        kw.update(keywords)
        if kw['defaults']:
            return kw  # pragma: no cover

        if fileformat == 'scorer':
            util.write_text_file(
                kw['filename'] + '.scorer', scorer2str(kw.get('scorer', self.rscorer)))
        else:
            self._output(fileformat, **kw)
