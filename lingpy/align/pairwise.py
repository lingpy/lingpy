# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-05 17:50
# modified : 2013-03-05 17:50
"""
Module provides classes and functions for pairwise alignment analyses.
"""

__author__="Johann-Mattis List"
__date__="2013-03-05"

# modules
from ..data import *
from ..sequence.sound_classes import *
try:
    from ..algorithm.cython import malign
    from ..algorithm.cython import calign
    from ..algorithm.cython import talign
except:
    from ..algorithm.cython import _malign as malign
    from ..algorithm.cython import _calign as calign
    from ..algorithm.cython import _talign as talign

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

    def __init__(
            self,
            seqs,
            seqB = False,
            **keywords
            ):

        # check, whether there are only two sequences or multiple sequence
        # pairs as input
        if seqB:
            self.seqs = [(seqs,seqB)]
        else:
            self.seqs = seqs

        # add the basic representation of sequences
        self.tokens = []
        self.prostrings = []

        # define a tokenizer function for convenience
        defaults = {
                "diacritics" : None,
                "vowels":None,
                "tones":None,
                "combiners":'\u0361\u035c',
                "breaks":'.-',
                "stress":"ˈˌ'",
                "merge_vowels" : True
                }
        for k in keywords:
            if k in defaults:
                defaults[k] = keywords[k]

        tokenize = lambda x: ipa2tokens(x,**defaults)

        # start to loop over data and create the stuff
        for k,(seqA,seqB) in enumerate(self.seqs):

            # get the tokens
            tokA,tokB = tokenize(seqA),tokenize(seqB)

            # get the prostrings
            proA,proB = prosodic_string(tokA),prosodic_string(tokB)

            # append the stuff
            self.tokens += [[tokA,tokB]]
            self.prostrings += [[proA,proB]]
    
    def __str__(self):

        try:
            a,b,c = self.alignments[0]
            out = '{0}\n{1}\n{2}'.format(
                    '\t'.join(a),
                    '\t'.join(b),
                    c
                    )
            for a,b,c in self.alignments[1:]:
                out += '\n\n' + '{0}\n{1}\n{2}'.format(
                        '\t'.join(a),
                        '\t'.join(b),
                        c
                        )
            return out

        # return tokens, if alignments aren't defined
        except:
            a,b = self.tokens[0]
            out = '{0}\n{1}'.format(
                    ''.join(a),
                    ''.join(b)
                    )
            for a,b in self.tokens[1:]:
                out += '\n\n' + '{0}\n{1}'.format(
                        ''.join(a),
                        ''.join(b)
                        )
            return out

    def __repr__(self):

        return str(self.seqs)
    
    def __len__(self):
        return len(self.seqs)
    
    def __getitem__(
            self,
            idx
            ):
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

    def _set_model(
            self,
            model = None
            ):
        """
        Define the sequence model for the calculation.

        Parameters
        ----------
        model : { None, lingpy.data.model.Model } (default=None)
            Specify the sound-class model to which the strings shall be
            converted.
        """
        
        if not model:
            self.model = sca
        else:
            self.model = model
        
        self.classes = []
        for clA,clB in map(
                lambda x:(tokens2class(x[0],self.model),tokens2class(x[1],self.model)),
                self.tokens
                ):
            self.classes += [(clA,clB)]

        self.weights = []
        for prA,prB in self.prostrings:
            self.weights += [(prosodic_weights(prA),prosodic_weights(prB))]
        
        self.scoredict = self.model.scorer

    def align(
            self,
            gop = -1,
            scale = 0.5,
            mode = 'global',
            factor = 0.3,
            restricted_chars = 'T_',
            distance = False,
            model = None,
            pprint = False
            ):
        """
        Align a pair of sequences or multiple sequence pairs.

        Parameters
        ----------
        gop : int (default=-1)
            The gap opening penalty (GOP).
        scale : float (default=0.5)
            The gap extension penalty (GEP), calculated with help of a scaling
            factor.
        mode : {'global','local','overlap','dialign'}
            The alignment mode, see :evobib:`List2012a` for details.
        factor : float (default = 0.3)
            The factor by which matches in identical prosodic position are
            increased.
        restricted_chars : str (default="T_")
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

        if hasattr(self,'model'):
            if not model:
                pass
            elif model == self.model:
                pass
            else:
                self._set_model(model)
        else:
            self._set_model(model)
        
        # create the alignments array
        self._alignments = []
        self.alignments = []
        
        # define local
        if mode == 'local':
            local = True
        else:
            local = False
        
        if not distance:
            self._alignments = calign.align_pairs(
                    self.classes,
                    self.weights,
                    self.prostrings,
                    gop,
                    scale,
                    factor,
                    self.scoredict,
                    mode,
                    restricted_chars
                    )
        else:
            self._alignments = calign.align_pairs(
                    self.classes,
                    self.weights,
                    self.prostrings,
                    gop,
                    scale,
                    factor,
                    self.scoredict,
                    mode,
                    restricted_chars,
                    distance = 1
                    )
        
        # switch back to alignments
        self.alignments = []
        for i,(almA,almB,sim) in enumerate(self._alignments):
            
            if mode != "local":
                self.alignments.append(
                        (
                            class2tokens(
                                self.tokens[i][0],almA
                                ),
                            class2tokens(
                                self.tokens[i][1],almB
                                ),
                            sim
                            )
                        )
            else:
                self.alignments.append(
                        (
                            class2tokens(
                                self.tokens[i][0],almA,local=True
                                ),
                            class2tokens(
                                self.tokens[i][1],almB,local=True
                                ),
                            sim
                            )
                        )

        # print the alignments, if this is chosen
        if pprint:
            print(self)
        else:
            pass

# the following functions provide solutions for convenience
def pw_align(
        seqA,
        seqB,
        gop = -1,
        scale = 0.5,
        scorer = False,
        mode = 'global',
        distance = False
        ):
    """
    Align two sequences in various ways.
    """

    # check whether the sequences are lists
    if type(seqA) == str or type(seqA) == tuple:
        seqA = list(seqA)
        seqB = list(seqB)
    elif type(seqA) != list:
        raise ValueError(
            "[!] Input sequences should be tuples, lists, or strings!"
            )
    if distance:
        distance = 1
    else:
        distance = 0
    
    if not scorer:
        scorer = {}
        for a in seqA:
            for b in seqB:
                if a == b:
                    scorer[a,b] = 1.0
                else:
                    scorer[a,b] = -1.0

    # start alignment
    return talign.align_pair(
            seqA,
            seqB,
            gop,
            scale,
            scorer,
            mode,
            distance
            )

def nw_align(
        seqA,
        seqB,
        scorer = False,
        gap = -1
        ):
    """
    Carry out the traditional Needleman-Wunsch algorithm.
    """
    # check whether the sequences are tuples
    if type(seqA) == str or type(seqA) == tuple:
        seqA = list(seqA)
        seqB = list(seqB)
    elif type(seqA) != list:
        raise ValueError(
            "[!] Input sequences should be tuples, lists, or strings!"
            )
    if not scorer:
        scorer = {}
        for a in seqA:
            for b in seqB:
                if a == b:
                    scorer[a,b] = 1.0
                else:
                    scorer[a,b] = -1.0

    return malign.nw_align(seqA,seqB,scorer,gap)

def edit_dist(
        seqA,
        seqB,
        normalized = False,
        restriction = ''
        ):
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
    restriction : {'cv'} (default='')
        Specify the restrictions to be used. Currently, only ``cv`` is
        supported. This prohibits matches of vowels with consonants.

    """
    # check whether the sequences are tuples
    if type(seqA) == str or type(seqA) == tuple:
        seqA = list(seqA)
        seqB = list(seqB)
    elif type(seqA) != list:
        raise ValueError(
            "[!] Input sequences should be tuples, lists, or strings!"
            )
    
    if restriction in ['cv','consonant-vowel']:
        
        resA = prosodic_string(seqA,'cv')
        resB = prosodic_string(seqB,'cv')

        return malign.restricted_edit_dist(
                seqA,
                seqB,
                resA,
                resB,
                normalized
                )
    
    return malign.edit_dist(seqA,seqB,normalized)

def sw_align(
        seqA,
        seqB,
        scorer = False,
        gap = -1
        ):
    """
    Carry out the traditional Smith-Waterman algorithm.
    """

    # check whether the sequences are tuples
    if type(seqA) == str or type(seqA) == tuple:
        seqA = list(seqA)
        seqB = list(seqB)
    elif type(seqA) != list:
        raise ValueError(
            "[!] Input sequences should be tuples, lists, or strings!"
            )
    if not scorer:
        scorer = {}
        for a in seqA:
            for b in seqB:
                if a == b:
                    scorer[a,b] = 1.0
                else:
                    scorer[a,b] = -1.0

    return malign.sw_align(seqA,seqB,scorer,gap)


def we_align(
        seqA,
        seqB,
        scorer = False,
        gap = -1
        ):
    """
    Carry out the traditional Waterman-Eggert algorithm.
    """

    # check whether the sequences are tuples
    if type(seqA) == str or type(seqA) == tuple:
        seqA = list(seqA)
        seqB = list(seqB)
    elif type(seqA) != list:
        raise ValueError(
            "[!] Input sequences should be tuples, lists, or strings!"
            )

    if not scorer:
        scorer = {}
        for a in seqA:
            for b in seqB:
                if a == b:
                    scorer[a,b] = 1.0
                else:
                    scorer[a,b] = -1.0

    return malign.we_align(seqA,seqB,scorer,gap)

def structalign(
        seqA,
        seqB
        ):
    """
    
    """
    return malign.structalign(seqA,seqB)

def turchin(
        seqA,
        seqB
        ):
    """
    Return cognate judgment based on Turchin et al.'s (2010) method.
    """
    if type(seqA) == str:
        seqA = ipa2tokens(seqA)
        seqB = ipa2tokens(seqB)
        
    classA = tokens2class(seqA,dolgo)
    classB = tokens2class(seqB,dolgo)

    if classA[0] == 'V': classA[0] = 'H'
    if classB[0] == 'V': classB[0] = 'H'

    if ''.join(classA).replace('V','')[:2] == ''.join(classB).replace('V','')[:2]:
        return 0
    else:
        return 1

