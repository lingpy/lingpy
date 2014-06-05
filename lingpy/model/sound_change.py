# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-05-11 12:16
# modified : 2014-05-11 12:16
"""
Module provides functions for modeling sound change.
"""

__author__="Johann-Mattis List"
__date__="2014-05-11"

from ..algorithm import misc
from ..sequence.sound_classes import bigrams, trigrams, pgrams, get_all_ngrams, class2tokens, ipa2tokens, prosodic_string
from ..align.pairwise import nw_align, Pairwise, sw_align, pw_align


class SoundChanger(object):
    
    def __init__(self, training, mode='pgrams', **keywords):
        """
        Parameters
        ----------
        training : list
            A list of tuples, containign the ancestral and the descendant
            words.
        """

        self.training = training
        self.alignments = Pairwise(training)(**keywords)
        self.ngrams = {}
        self.converter = {}

        self.sources = [t[0] for t in training]
        self.targets = [t[1] for t in training]
        
        max_ngram = 0

        for almA,almB,sim in self.alignments:
            
            if mode == 'pgrams':
                # we need to detect gaps in target first, otherwise we cannot
                # determine the "correct" prosodic strings in the beginning,
                # since gaps are strangely encoded in prostrings
                gaps = [i for i in range(len(almA)) if almA[i] == '-']

                # convert alignments to programs
                pA = prosodic_string([almA[i] for i in range(len(almA)) if i
                    not in gaps], **keywords)
                pA = class2tokens(pA,almA)
                progA = list(zip(almA,pA))
                #progA = pgrams([almA[i] for i in range(len(almA)) if i not in
                #    gaps],**keywords)

                progB = pgrams(almB,**keywords)

            elif mode == 'simple':
                progA = list(zip(almA,['#']+almA[1:-1]+['$']))
                progB = list(zip(almB,['#']+almB[1:-1]+['$']))

            elif mode == 'bigrams':
                progA = bigrams(almA)
                progB = bigrams(almB)

            
            # zip source and target
            zips = list(zip(progA,progB))

            # compute all ngrams
            ngrams = get_all_ngrams(zips)

            # iterate over all ngrams and add them to the mother dictionary
            for ngram in ngrams:

                # get upper and lower
                upper = tuple([g[0] for g in ngram if g[0][0] != '-'])
                lower = tuple([g[1] for g in ngram if g[1][0] != '-'])

                # backup for pgrams
                #if mode == 'pgrams':
                #    upper = tuple(pgrams(ipa2tokens([g[0] for g in upper]),
                #        **keywords))

                # append to ngrams dictionary
                try:
                    self.ngrams[upper,lower] += 1
                except:
                    self.ngrams[upper,lower] = 1

                try:
                    self.converter[upper] += [lower]
                except:
                    self.converter[upper] = [lower]

                if len(upper) > max_ngram:
                    max_ngram = len(upper)
        
        # carry out post-processing of ngrams, delete all cases where ngrams
        # occur only once
        wunnies = [ngram for ngram in self.ngrams if self.ngrams[ngram] == 1
                and len(self.converter[ngram[0]]) == 1]
        
        full_words = []
        part_words = []
        for wun in wunnies:
            word = ''.join([w[0] for w in wun[0]])
            if word in self.sources:
                full_words += [word]
            else:
                part_words += [(wun,word)]
        for wun,word in part_words:
            tmp = [w for w in full_words if word in w]
            if tmp:
                del self.ngrams[wun]

        self._max_ngram = max_ngram

    def get_trivials(self):
        """
        Get trivial sound change pairs.
        """
        self.trivials = []
        
        for i in range(1,self._max_ngram+1):
            
            sources = [k[0] for k in self.ngrams if len(k[0]) == i]

            for j,source in enumerate(sources):

                if sources.count(source) <= 1:

                    if self._split_string(source):
                        pass
                    else:
                        self.trivials += [source]
        
    def _split_string(self, sequence, debug=False):
        """
        Dynamic-programming approach to splitting all input strings into substrings which have unique counterparts in the target language.
        """
        queue = [[tuple([]),-1]]

        while queue:

            current,idx = queue.pop(0)
            if debug:
                print(idx, current, current in self.trivials)

            if idx == len(sequence)-1:
                if current in self.trivials:
                    return True
            
            else:
                subseq = current+tuple([sequence[idx+1]])
                
                if subseq in self.trivials:
                    queue += [[tuple([]),idx+1]]
                    queue += [[subseq,idx+1]]
                else:
                    queue += [[subseq, idx+1]]
        
        return False

    
    def transform(self, sequence, mode='pgrams', debug=False, **keywords):
        
        if mode == 'pgrams':
            seq = tuple(pgrams(sequence, **keywords))

        elif mode == 'simple':
            seq = ipa2tokens(sequence)
            seq = tuple(zip(seq,['#']+seq[1:-1]+['$']))

        elif mode == 'bigrams':
            seq = bigrams(ipa2tokens(sequence))

        queue = [[tuple([]),-1,tuple([])]]
        output = []

        while queue:

            current,idx,outseq = queue.pop(0)

            if debug:
                print(idx,current,outseq,current in self.trivials)

            if idx == len(seq) -1:
                if current in self.trivials:
                    output += [outseq+self.converter[current][0]]
                else:
                    pass
            else:
                subseq = current + tuple([seq[idx+1]])

                if subseq in self.trivials:
                    queue += [[tuple([]),idx+1,outseq+self.converter[subseq][0]]]
                    queue += [[subseq,idx+1,outseq]]
                else:
                    queue += [[subseq, idx+1, outseq]]

        routput = []
        for p in output:
            w = ''.join([k[0] for k in p])
            if w not in routput:
                routput += [w]
        return routput

        
