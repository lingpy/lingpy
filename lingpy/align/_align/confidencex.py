# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-02-18 12:26
# modified : 2014-02-18 12:26
"""
Calculate confidence scores for the scoring functions in alignment plots.
"""

__author__="Johann-Mattis List"
__date__="2014-02-18"

from ...sequence.sound_classes import class2tokens

def get_confidence(alms, scorer, ref='lexstatid', gap_weight=1):
    """
    Function creates confidence scores for a given set of alignments.

    Parameters
    ----------
    alms : :py:class`~lingpy.align.sca.Alignments`
        An *Alignments* object containing already aligned strings.
    scorer : :py:class:`~lingpy.algorithm._misc.ScoreDict`
        A *ScoreDict* object which gives similarity scores for all segments in
        the alignment.
    ref : str (default="lexstatid")
        The reference entry-type, referring to the cognate-set to be used for
        the analysis.
    """
    
    # store all values for average scores
    values = []

    for key,msa in alms.msa[ref].items():
        
        # get basic stuff
        idxs = msa['ID']
        taxa = msa['taxa']

        # get numerical representation of alignments
        if scorer:
            alignment = [class2tokens(
                alms[idxs[i], 'numbers'],
                msa['alignment'][i]) for i in range(len(idxs))]
        else:
            alignment = msa['alignment']
        
        # create new array for confidence
        confidence_matrix = []

        # iterate over each taxon
        for i,taxon in enumerate(taxa):
            idx = alms.taxa.index(taxon)+1

            # get the numerical sequence
            nums = alignment[i]
            
            # store confidences per line
            confidences = []

            # iterate over the sequence
            for j,num in enumerate(nums):

                col = [alm[j] for alm in alignment]
                
                score = 0
                count = 0
                    
                for k,numB in enumerate(col):
                    if k != i:
                        if num == '-' and numB == '-':
                            pass
                        else:
                            gaps = False
                            if num == '-' and numB != '-':
                                numA = str(idx)+'.X.-'
                                gaps = True
                            elif numB == '-' and num != '-':
                                numB = str(alms.taxa.index(taxa[k]))+'.X.-'
                                numA = num
                                gaps = True
                            else:
                                numA = num
                            
                            scoreA = scorer[numA,numB]
                            scoreB = scorer[numB,numA]
                            this_score = max(scoreA,scoreB)
                            
                            if not gaps:
                                score += this_score
                                count += 1
                            else:
                                score += this_score * gap_weight
                                count += gap_weight
                            
                if count:
                    score = score / count #(len(col) - gaps * gap_weight)
                else:
                    score = -25

                confidences += [int(score+0.5)]
                values += [int(score+0.5)]
            confidence_matrix += [confidences]

        # append confidence matrix to alignments
        alms.msa[ref][key]['confidence'] = confidence_matrix


    # sort the values
    values = sorted(set(values))

    # make conversion to scale of 100 values
    converter = {}
    step = 100 / (len(values)+1)
    for i,score in enumerate(values):
        converter[score] = int(step * score+0.5)

    # iterate over keys again
    for key,msa in alms.msa[ref].items():
        
        # get basic stuff
        for i,line in enumerate(msa['confidence']):
            for j,cell in enumerate(line):
                msa['confidence'] = converter[cell]



