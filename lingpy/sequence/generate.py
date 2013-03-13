# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-13 08:02
# modified : 2013-03-13 08:02
"""
Module provides simple basic classes for sequence generation using Markov models.
"""

__author__="Johann-Mattis List"
__date__="2013-03-13"

# builtin
import random

# internal
from ..data import *
from .sound_classes import *

class MCBasic(object):
    """
    Basic class for creating Markov chains from sequence training data.

    Parameters
    ----------
    seq : list
        A list of sequences. Sequences are assumed to be tokenized, i.e. they
        should be either passed as lists or as tuples.

    """

    def __init__(
            self,
            seqs
            ):

        self.seqs = seqs
        
        # create distribution
        self.dist = {}
        for seq in self.seqs:

            for s1,s2 in zip(['#']+seq,seq+['$']):
                try:
                    self.dist[s1] += [s2]
                except:
                    self.dist[s1] = [s2]

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

class MCPhon(MCBasic):
    """
    Class for the creation of phonetic sequences ("pseudo words").

    Parameters
    ----------
    words : list
        List of phonetic sequences. This list can contain tokenized
        sequences (lists or tuples), or simple untokenized IPA strings.

    tokens : bool (default=False)
        If set to True, no tokenization of input sequences is carried out.

    prostring : list (default=[])
        List containing the prosodic profiles of the input sequences. If the
        list is empty, the profiles are generated automatically.

    """
    
    def __init__(
            self,
            words,
            tokens=False,
            prostrings=[]
            ):
        
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
        Generate a string from the Markov chain created from the training data.
        
        Parameters
        ----------
        new : bool (default=True)
            Determine whether the string created should be different from the
            training data or not.
        tokens : bool (default=False)
            If set to *True* he full list of tokens that was internally used to
            represent the sequences as a Markov chain is returned.
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

