"""
Basic module for pairwise alignment analyses.
"""
# modules
from __future__ import division,print_function
from ..data import *
from ..algorithm import *


class _Pairwise(object):
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

    merge_vowels : bool (default=True)
        Indicate, whether neighboring vowels should be merged into
        diphtongs, or whether they should be kept separated during the
        analysis.

    """
    def __init__(
            self,
            seqs,
            seqB = False,
            merge_vowels = True
            ):

        # check, whether there are only two sequences or multiple sequence
        # pairs as input
        if seqB:
            self.seqs = [(seqs,seqB)]
        else:
            self.seqs = seqs
            
        # create a tokenized representation of all sequences
        self.tokens = []
        self.sonars = []
        self.numbers = []
        self.prosodics = []
        for k,(i,j) in enumerate(self.seqs):
            seqA = ipa2tokens(i,merge_vowels=merge_vowels)
            seqB = ipa2tokens(j,merge_vowels=merge_vowels)
            sonA = [int(x) for x in tokens2class(seqA,art)]
            sonB = [int(x) for x in tokens2class(seqB,art)]
            self.tokens.append([seqA,seqB])
            self.sonars.append([sonA,sonB])
            self.numbers.append([
                [str(k) +'.0.'+ str(p) for p in range(1,len(seqA)+1)],
                [str(k) +'.1.'+ str(p) for p in range(1,len(seqB)+1)]
                ])
            self.prosodics.append([
                prosodic_string(sonA),
                prosodic_string(sonB)
                ])
    def __str__(self):

        # check for alignments
        #@todo: check bug in output!
        try:
            a,b,c = self.alignments[0]
            out = '{0}\n{1}\n{2}'.format(
                    '\t'.join(a).encode('utf-8'),
                    '\t'.join(b).encode('utf-8'),
                    c)
            for a,b,c in self.alignments[1:]:
                out += '\n\n' + '{0}\n{1}\n{2}'.format(
                        '\t'.join(a).encode('utf-8'),
                        '\t'.join(b).encode('utf-8'),
                        c)
            return out
        # return tokens, if alignments aren't defined
        except:
            a,b = self.tokens[0]
            out = '{0}\n{1}'.format(
                    '\t'.join(a).encode('utf-8'),
                    '\t'.join(b).encode('utf-8'))
            for a,b in self.tokens[1:]:
                out += '\n\n' + '{0}\n{1}'.format(
                        '\t'.join(a).encode('utf-8'),
                        '\t'.join(b).encode('utf-8'))
            return out
    
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
            model = 'sca'
            ):
        """
        Define the sequence model for the calculation.
        """
        
        try:
            self.model = eval(model)
        except:
            self.model = model
        
        self.classes = [[tokenA[:],tokenB[:]] for tokenA,tokenB in self.tokens]
        
        letterA,letterB = [],[]
        
        # mind that the C++ function expects lists as input.
        for i,(tokenA,tokenB) in enumerate(self.classes):
            self.classes[i][0] = list(tokens2class(tokenA,self.model))
            self.classes[i][1] = list(tokens2class(tokenB,self.model))
            letterA += self.classes[i][0]
            letterB += self.classes[i][1]

        self.scoredict = self.model.scorer

    def _get(self,number,value='tokens',error=(0,'-')):
        """
        Return a specific segment of the input sequences.
        """
        
        if number == error[1]:
            return error[1]

        number = [int(i) for i in number.split('.')]

        if value == 'tokens':
            return self.tokens[number[0]][number[1]][number[2]-1]
        elif value == 'classes':
            return self.classes[number[0]][number[1]][number[2]-1]

    def _align_pairwise(
            self,
            mode = 'global',
            scale = (1.2,1.0,1.1),
            factor = 0.3,
            gop = -2,
            gep_scale = 0.5,
            restricted_chars = '_T'
            ):

        alms = []
        weights = []
        for seqA,seqB in self.prosodics:
            prosA = gop * array(
                    prosodic_weights(
                        seqA,
                        scale,
                        factor
                        ),
                    dtype='float'
                    )
            prosB = gop * array(
                    prosodic_weights(
                        seqB,
                        scale,
                        factor
                        ),
                    dtype = 'float'
                    )
            weights.append([prosA.tolist(),prosB.tolist()])
        restrictions = []
        prosodic_strings = []
        for prosA,prosB in self.prosodics:
            tmpA,tmpB = [],[]
            for i,char in enumerate(prosA):
                if char in restricted_chars:
                    tmpA.append(-(i+1))
                else:
                    tmpA.append(i+1)
            for i,char in enumerate(prosB):
                if char in restricted_chars:
                    tmpB.append(-(i+1))
                else:
                    tmpB.append(i+1)

            restrictions.append([tmpA,tmpB])
            prosodic_strings.append([prosA,prosB])
        
        self._aligned_classes = align_sequence_pairs(
                self.classes,
                weights,
                restrictions,
                prosodic_strings,
                self.scoredict,
                gep_scale,
                0.25 * factor,
                mode
                )
        self.weights = weights

        # the following lines convert the output to the different formats
        self._aligned_numbers = []
        self.alignments = []
        for i,(a,b,c) in enumerate(self._aligned_classes):
            numA = []
            numB = []
            k = 0
            for j in range(len(a)):
                if a[j] not in '-*':
                    numA.append(str(i)+'.0.'+str(k+1))
                    k += 1
                elif a[j] == '-':
                    numA.append('-')
                else:
                    k += 1

            k = 0
            for j in range(len(b)):
                if b[j] not in '-*':
                    numB.append(str(i)+'.1.'+str(k+1))
                    k += 1
                elif b[j] == '-':
                    numB.append('-')
                else:
                    k += 1
            self._aligned_numbers.append((numA,numB,c))
            self.alignments.append(
                    [
                        [self._get(k) for k in numA],
                        [self._get(k) for k in numB],
                        c
                        ]
                    )

    def align(
            self,
            model = 'sca',
            mode = 'global',
            gop = -5,
            gep_scale = float(0.5),
            scale = None, #(float(1.2),float(1.0),float(1.1)),
            factor = float(0.3),
            restricted_chars = '_T'
            ):
        """
        Align two sequences or a list of sequence pairs pairwise.

        Parameters
        ----------

        model : string (default="sca")
            A string indicating the name of the
            :py:class:`~lingpy.data.model.Model` object that shall be used for
            the analysis.
            Currently, three models are supported:
            
            * "dolgo" -- a sound-class model based on :evobib:`Dolgopolsky1986`,

            * "sca" -- an extension of the "dolgo" sound-class model based on
              :evobib:`List2012a`, and
            
            * "asjp" -- an independent sound-class model which is based on the
              sound-class model of :evobib:`Brown2008` and the empirical data of
              :evobib:`Brown2011`.

        mode : string (default="global")
            A string indicating which kind of alignment analysis should be
            carried out. Select between: 
            
            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1971`,

            * "local" -- local alignment analysis based on the Smith-Waterman
              algorithm :evobib:`Smith1981`,

            * "overlap" -- overlap alignment analysis where gaps introduced in
              the beginning or the end of the alignment are scored with 0.

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.
        
        gop : int (default=-5)
            The gap opening penalty (gop) on which the analysis shall be based.

        gep_scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1982`.

        scale : tuple or list (default=(3,1,2))
            The scaling factors for the modificaton of gap weights. The first
            value corresponds to sites of ascending sonority, the second value
            to sites of maximum sonority, and the third value corresponds to
            sites of decreasing sonority.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.

        restricted_chars : string (default="T")
            Define which characters of the prosodic string of a sequence
            reflect its secondary structure (cf. :evobib:`List2012a`) and
            should therefore be aligned specifically. This defaults to "T",
            since this is the character that represents tones in the prosodic
            strings of sequences.
              
        """

        self._set_model(model)
        self._align_pairwise(
                mode,
                scale,
                factor,
                gop,
                gep_scale,
                restricted_chars
                )
    
    def _self_score(
            self,
            seq
            ):
        """
        Return the score of a sequence with itself.
        """

        s = 0
        for x in seq:
            if x != '-':
                s += self.model.scorer[x,x]
        return s

    def distance(
            self
            ):
        """
        Return the distance from the alignments.
        """
        
        for i,(almA,almB,sAB) in enumerate(self._aligned_classes):
            sA,sB = self._self_score(almA),self._self_score(almB)
            d = 1 - ((2 * sAB) / (sA+sB))
            self.alignments[i][2] = d
