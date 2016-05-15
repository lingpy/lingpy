# *-* coding: utf-8 *-*
"""
Module provides classes and functions for pairwise alignment analyses.
"""
from __future__ import division, print_function, unicode_literals

from six import text_type
from lingpy.util import setdefaults, multicombinations2, as_string 
from lingpy.settings import rcParams
from lingpy.sequence.sound_classes import (
    ipa2tokens, prosodic_string, tokens2class, prosodic_weights, class2tokens,
)
from lingpy.algorithm import malign
from lingpy.algorithm import calign
from lingpy.algorithm import talign


class Pairwise(object):
    """
    Basic class for the handling of pairwise sequence alignments (PSA).

    Parameters
    ----------

    seqs : string list
        Either the first string of a sequence pair that shall be aligned,
        or a list of sequence tuples.

    seqB : string or bool (default=None)
        Define the second sequence that shall be aligned with the first
        sequence, if only two sequences shall be compared.

    """

    def __init__(self, seqs, seqB=False, **keywords):
        # check, whether there are only two sequences or multiple sequence
        # pairs as input
        if seqB:
            self.seqs = [(seqs, seqB)]
        else:
            self.seqs = seqs

        # add the basic representation of sequences
        self.tokens = []
        self.prostrings = []

        # define a tokenizer function for convenience
        defaults = {
            "diacritics": rcParams['diacritics'],
            "vowels": rcParams['vowels'],
            "tones": rcParams['tones'],
            "combiners": rcParams['combiners'],
            "breaks": rcParams['breaks'],
            "stress": rcParams['stress'],
            "merge_vowels": rcParams['merge_vowels']
        }
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        tokenize = lambda x: ipa2tokens(x, **keywords)

        # start to loop over data and create the stuff
        for k, (seqA, seqB) in enumerate(self.seqs):
            # get the tokens
            tokA, tokB = tokenize(seqA), tokenize(seqB)

            # get the prostrings
            proA, proB = \
                prosodic_string(tokA, **keywords), prosodic_string(tokB, **keywords)

            # append the stuff
            self.tokens += [[tokA, tokB]]
            self.prostrings += [[proA, proB]]

    def __str__(self):
        try:
            a, b, c = self.alignments[0]
            out = '{0}\n{1}\n{2}'.format('\t'.join(a), '\t'.join(b), c)
            for a, b, c in self.alignments[1:]:
                out += '\n\n' + '{0}\n{1}\n{2}'.format('\t'.join(a), '\t'.join(b), c)
            return out

        # return tokens, if alignments aren't defined
        except:
            a, b = self.tokens[0]
            out = '{0}\n{1}'.format(''.join(a), ''.join(b))
            for a, b in self.tokens[1:]:
                out += '\n\n' + '{0}\n{1}'.format(''.join(a), ''.join(b))
            return out

    def __call__(self, **keywords):
        self.align(**keywords)
        return self.alignments

    def __repr__(self):
        return text_type(self.seqs)

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        """
        Return specified values.
        """
        try:
            data = idx[1]
            idx = idx[0]
        except:
            data = 'w'

        if data == 'w':
            return self.seqs[idx]
        elif data == 'c':
            return self.classes[idx]
        elif data == 't':
            return self.tokens[idx]
        elif data == 'a':
            return self.alignments[idx]

    def _set_model(self, **keywords):
        """
        Define the sequence model for the calculation.

        Parameters
        ----------
        model : { None, Model } (default=None)
            Specify the sound-class model to which the strings shall be
            converted.
        """
        defaults = dict(
            model=rcParams['sca'],
            stress=rcParams['stress'],
            transform=rcParams['align_transform'])
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        if isinstance(keywords['model'], (text_type, str)):
            self.model = rcParams[keywords['model']]
        else:
            self.model = keywords['model']

        self.classes = []
        for clA, clB in map(
            lambda x: (
                tokens2class(x[0], self.model, stress=keywords['stress']),
                tokens2class(x[1], self.model, stress=keywords['stress'])),
            self.tokens
        ):
            self.classes += [(clA, clB)]

        self.weights = []
        for prA, prB in self.prostrings:
            self.weights += [(
                prosodic_weights(prA, _transform=keywords['transform']),
                prosodic_weights(prB, _transform=keywords['transform'])
            )]

        self.scoredict = self.model.scorer

    def align(self, **keywords):
        """
        Align a pair of sequences or multiple sequence pairs.

        Parameters
        ----------
        gop : int (default=-1)
            The gap opening penalty (GOP).
        scale : float (default=0.5)
            The gap extension penalty (GEP), calculated with help of a scaling
            factor.
        mode : {"global","local","overlap","dialign"}
            The alignment mode, see :evobib:`List2012a` for details.
        factor : float (default = 0.3)
            The factor by which matches in identical prosodic position are
            increased.
        restricted_chars : str (default="T\_")
            The restricted chars that function as an indicator of syllable or
            morpheme breaks for secondary alignment, see :evobib:`List2012c`
            for details.
        distance : bool (default=False)
            If set to *True*, return the distance instead of the similarity
            score. Distance is calculated using the formula by
            :evobib:`Downey2008`.
        model : { None, ~lingpy.data.model.Model }
            Specify the sound class model that shall be used for the analysis.
            If no model is specified, the default model of :evobib:`List2012a`
            will be used.
        pprint : bool (default=False)
            If set to *True*, the alignments are printed to the screen.

        """
        setdefaults(
            keywords,
            gop=-1,
            scale=0.5,
            mode='global',
            factor=0.3,
            restricted_chars='T_',
            distance=False,
            model=rcParams['sca'],
            pprint=False,
            transform=rcParams['align_transform'])

        if hasattr(self, 'model'):
            if keywords['model'] != self.model:
                self._set_model(**keywords)
        else:
            self._set_model(**keywords)

        # create the alignments array
        self._alignments = calign.align_pairs(
            self.classes,
            self.weights,
            self.prostrings,
            keywords['gop'],
            keywords['scale'],
            keywords['factor'],
            self.scoredict,
            keywords['mode'],
            keywords['restricted_chars'],
            distance=1 if keywords['distance'] else 0)

        # switch back to alignments
        self.alignments = []
        for i, (almA, almB, sim) in enumerate(self._alignments):
            self.alignments.append((
                class2tokens(self.tokens[i][0], almA, local=keywords['mode'] == "local"),
                class2tokens(self.tokens[i][1], almB, local=keywords['mode'] == "local"),
                sim))

        # print the alignments, if this is chosen
        as_string(self, pprint=keywords['pprint'])

def _get_scorer(seqA, seqB):
    scorer = {}
    for a in seqA:
        for b in seqB:
            scorer[a, b] = 1.0 if a == b else -1.0
    return scorer


# the following functions provide solutions for convenience
def pw_align(
        seqA,
        seqB,
        gop=-1,
        scale=0.5,
        scorer=False,
        mode='global',
        distance=False,
        **keywords):
    """
    Align two sequences in various ways.

    Parameters
    ----------
    seqA, seqB : {text_type, list, tuple}
        The input strings. These should be iterables, so you can use tuples,
        lists, or strings.
    scorer : dict (default=False)
        If set to c{False} a scorer will automatically be calculated,
        otherwise, the scorer needs to be passed as a dictionary that covers
        all segment matches between the input strings.
    gop : int (default=-1)
        The gap opening penalty.
    scale : float (default=0.5)
        The gap extension scale. This scale is similar to the gap extension
        penalty, but in contrast to the traditional GEP, it "scales" the gap
        opening penalty.
    mode : {"global", "local", "diagonal", "overlap"} (default="global")
        Select between one of the four different alignment modes regularly
        implemented in LingPy, see :evobib:`List2012a` for details.
    distance : bool (default=False)
        If set to c{True} return the distance score following the formula by
        :evobib:`Downey2008`. Otherwise, return the basic similarity score.

    Examples
    --------
    Align two words using the dialign algorithm::
        >>> seqA = 'fat cat'
        >>> setB = 'catfat'
        >>> pw_align(seqA, seqB, mode='dialign')
        (['f', 'a', 't', ' ', 'c', 'a', 't', '-', '-', '-'],
         ['-', '-', '-', '-', 'c', 'a', 't', 'f', 'a', 't'],
         3.0)

    """

    # check whether the sequences are lists
    if isinstance(seqA, (text_type, tuple)):
        seqA = list(seqA)
        seqB = list(seqB)
    elif not isinstance(seqA, list):
        raise ValueError("Input should be tuple, list, or string.")

    distance = 1 if distance else 0

    if not scorer and distance == 0:
        scorer = _get_scorer(seqA, seqB)
    elif not scorer and distance == 1:
        scorer = {}
        for (i, a), (j, b) in multicombinations2(enumerate(sorted(set(seqA + seqB)))):
            scorer[a, b] = 1.0 if a == b else -1.0
            scorer[b, a] = 1.0 if a == b else -1.0

    # start alignment
    return talign.align_pair(seqA, seqB, gop, scale, scorer, mode, distance)


def nw_align(seqA, seqB, scorer=False, gap=-1):
    """
    Carry out the traditional Needleman-Wunsch algorithm.

    Parameters
    ----------
    seqA, seqB : {str, list, tuple}
        The input strings. These should be iterables, so you can use tuples,
        lists, or strings.
    scorer : dict (default=False)
        If set to c{False} a scorer will automatically be calculated,
        otherwise, the scorer needs to be passed as a dictionary that covers
        all segment matches between the input strings (segment matches need to
        be passed as tuples of two segments, following the order of the input
        sequences). Note also that the scorer can well be asymmetric, so you
        could also use it for two completely different alphabets. All you need
        to make sure is that the tuples representing the segment matches follow
        the order of your input sequences.
    gap : int (default=-1)
        The gap penalty.

    Notes
    -----
    The Needleman-Wunsch algorithm (see :evobib:`Needleman1970`) returns a global
    alignment of two sequences.

    Returns
    -------
    alm : tuple
        A tuple consisting of the aligments of the first and the second
        sequence, and the alignment score.

    Examples
    --------

    >>> seqA = 'fat cat'
    >>> setB = 'catfat'
    >>> nw_align(seqA,seqB)
    (['f', 'a', 't', ' ', 'c', 'a', 't'], ['c', 'a', 't', '-', 'f', 'a', 't'], 1)

    Use your own scorer (make sure all characters are covered, or you use a
    default dict). We start with a scorer that is "normal", with identical symbols getting
    identical scores:

    >>> scorer = { ('a','a'): 1, ('a','b'):-1, ('b','a'):-1, ('b', 'b'): 1}
    >>> seqA, seqB = 'abab', 'baba'
    >>> almA, almB, sim = nw_align(seqA, seqB, scorer=scorer)
    >>> print(' '.join(almA)+'\n'+' '.join(almB), "(sim={0})".format(sim))
    a b a b -
    - b a b a (sim=1)

    Nothing unexpected so far, you could reach the same result without the
    scorer. But now let's make a scorer that favors mismatches for our little
    two-letter alphabet.

    >>> scorer = { ('a','b'): 1, ('a','a'):-1, ('b','b'):-1, ('b', 'a'): 1}
    >>> seqA, seqB = 'abab', 'baba'
    >>> almA, almB, sim = nw_align(seqA, seqB, scorer=scorer)
    >>> print(' '.join(almA)+'\n'+' '.join(almB), "(sim={0})".format(sim))
    a b a b
    b a b a (sim=4)

    Now, let's analyse two strings which are completely different, but where we
    use the scorer to define mappings between the segments. We simply do this
    by using lower case letters in one and upper case letters in the other
    case, which will, of course, be treated as different symbols in Python:

    >>> scorer = { ('A','a'): 1, ('A','b'):-1, ('B','a'):-1, ('B', 'B'): 1}
    >>> seqA, seqB = 'ABAB', 'aa'
    >>> almA, almB, sim = nw_align(seqA, seqB, scorer=scorer)
    >>> print(' '.join(almA)+'\n'+' '.join(almB), "(sim={0})".format(sim))
    A B A B
    a - a - (sim=0)

    """
    # check whether the sequences are tuples
    if isinstance(seqA, (text_type, tuple)):
        seqA = list(seqA)
        seqB = list(seqB)
    elif not isinstance(seqA, list):
        raise ValueError("Input should be tuple, list, or string.")

    return malign.nw_align(seqA, seqB, scorer or _get_scorer(seqA, seqB), gap)


def edit_dist(seqA, seqB, normalized=False, restriction=''):
    """
    Return the edit distance between two strings.

    Parameters
    ----------
    seqA,seqB : str
        The strings that shall be compared.
    normalized : bool (default=False)
        Specify whether the normalized edit distance shall be returned. If no
        restrictions are chosen, the edit distance is normalized by dividing by
        the length of the longer string. If *restriction* is set to *cv*
        (consonant-vowel), the edit distance is normalized by the length of the
        alignment.
    restriction : {"cv"} (default="")
        Specify the restrictions to be used. Currently, only ``cv`` is
        supported. This prohibits matches of vowels with consonants.

    Notes
    -----
    The edit distance was first formally defined by V. I. Levenshtein
    (:evobib:`Levenshtein1965`). The first algorithm to compute the edit
    distance was proposed by Wagner and Fisher (:evobib:`Wagner1974`).

    Returns
    -------
    dist : {int float}
        The edit distance, which is a float if normalized is set to c{True},
        and an integer otherwise.

    Examples
    --------
    >>> seqA = 'fat cat'
    >>> setB = 'catfat'
    >>> edit_dist(seqA, seqB)
    3

    """
    # check whether the sequences are tuples
    if isinstance(seqA, (text_type, tuple)):
        seqA = list(seqA)
        seqB = list(seqB)
    elif not isinstance(seqA, list):
        raise ValueError("Input should be tuple, list, or string.")

    if restriction in ['cv', 'consonant-vowel']:
        resA = prosodic_string(seqA, 'cv')
        resB = prosodic_string(seqB, 'cv')
        return malign.restricted_edit_dist(seqA, seqB, resA, resB, normalized)

    return malign.edit_dist(seqA, seqB, normalized)


def sw_align(seqA, seqB, scorer=False, gap=-1):
    """
    Carry out the traditional Smith-Waterman algorithm.

    Parameters
    ----------
    seqA, seqB : {str, list, tuple}
        The input strings. These should be iterables, so you can use tuples,
        lists, or strings.
    scorer : dict (default=False)
        If set to c{False} a scorer will automatically be calculated,
        otherwise, the scorer needs to be passed as a dictionary that covers
        all segment matches between the input strings.
    gap : int (default=-1)
        The gap penalty.

    Notes
    -----
    The Smith-Waterman algorithm (see :evobib:`Smith1981`) returns a local
    alignment between two sequences. A local alignment is an alignment of those
    subsequences of the input sequences that yields the highest score.

    Returns
    -------
    alm : tuple
        A tuple consisting of prefix, alignment, and suffix of the first and
        the second sequence, and the alignment score.

    Examples
    --------

    >>> seqA = 'fat cat'
    >>> setB = 'catfat'
    >>> sw_align(seqA, seqB)
    (([], ['f', 'a', 't'], [' ', 'c', 'a', 't']),
     (['c', 'a', 't'], ['f', 'a', 't'], []),
     3.0)
    """

    # check whether the sequences are tuples
    if isinstance(seqA, (text_type, tuple)):
        seqA = list(seqA)
        seqB = list(seqB)
    elif not isinstance(seqA, list):
        raise ValueError("Input should be tuple, list, or string.")

    return malign.sw_align(seqA, seqB, scorer or _get_scorer(seqA, seqB), gap)


def we_align(seqA, seqB, scorer=False, gap=-1):
    """
    Carry out the traditional Waterman-Eggert algorithm.

    Parameters
    ----------
    seqA, seqB : {str, list, tuple}
        The input strings. These should be iterables, so you can use tuples,
        lists, or strings.
    scorer : dict (default=False)
        If set to c{False} a scorer will automatically be calculated,
        otherwise, the scorer needs to be passed as a dictionary that covers
        all segment matches between the input strings.
    gap : int (default=-1)
        The gap penalty.

    Notes
    -----
    The Waterman-Eggert algorithm (see :evobib:`Waterman1987`) returns *all*
    local matches between two sequences.

    Returns
    -------
    alms : list
        A list consisting of tuples. Each tuple gives the alignment of one of
        the subsequences of the input sequences. Each tuple contains the
        aligned part of the first, the aligned part of the second sequence, and
        the score of the alignment.

    Examples
    --------

    >>> seqA = 'fat cat'
    >>> setB = 'catfat'
    >>> we_align(seqA, seqB)
    [(['f', 'a', 't'], ['f', 'a', 't'], 3.0),
     (['c', 'a', 't'], ['c', 'a', 't'], 3.0)]

    """

    # check whether the sequences are tuples
    if isinstance(seqA, (text_type, tuple)):
        seqA = list(seqA)
        seqB = list(seqB)
    elif not isinstance(seqA, list):
        raise ValueError("Input should be tuple, list, or string.")

    return malign.we_align(seqA, seqB, scorer or _get_scorer(seqA, seqB), gap)


def structalign(seqA, seqB):
    """
    Experimental function for testing structural alignment algorithms.
    """
    return malign.structalign(seqA, seqB)


def turchin(seqA, seqB, model='dolgo', **keywords):
    """
    Return cognate judgment based on the method by :evobib:`Turchin2010`.

    Parameters
    ----------
    seqA, seqB : {str, list, tuple}
        The input strings. These should be iterables, so you can use tuples,
        lists, or strings.
    model : {"asjp", "sca", "dolgo"} (default="dolgo")
        A sound-class model instance or a string that denotes one of the
        standard sound class models used in LingPy.

    Returns
    -------
    cognacy : {0, 1}
        The cognacy assertion which is either 0 (words are probably cognate) or
        1 (words are not likely to be cognate).

    """
    if text_type(model) == model:
        model = rcParams[model]
    elif hasattr(model, 'info'):
        pass
    else:
        raise ValueError("[!] No valid model instance selected.")

    if isinstance(seqA, (text_type, str)):
        seqA = ipa2tokens(seqA)
        seqB = ipa2tokens(seqB)

    classA = tokens2class(seqA, model)
    classB = tokens2class(seqB, model)

    if classA[0] in model.vowels:
        classA[0] = 'H'
    if classB[0] in model.vowels:
        classB[0] = 'H'

    if ''.join([k for k in classA if k not in model.vowels])[:2] == \
            ''.join([k for k in classB if k not in model.vowels])[:2]:
        return 0
    else:
        return 1
