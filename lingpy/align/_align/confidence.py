from __future__ import division, unicode_literals
"""
Calculate confidence scores for the scoring functions in alignment plots.
"""
import cgi

from lingpy.sequence.sound_classes import class2tokens, token2class
from lingpy.settings import rcParams
from lingpy.util import charstring, dotjoin


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

    # store all correspondences
    corrs = {}

    # store occurrences
    occs = {}

    for key, msa in alms.msa[ref].items():
        # get basic stuff
        idxs = msa['ID']
        taxa = msa['taxa']
        concept = cgi.escape(alms[idxs[0], 'concept'], True)

        # get numerical representation of alignments
        if scorer:
            alignment = [class2tokens(
                alms[idxs[i], 'numbers'],
                msa['alignment'][i]) for i in range(len(idxs))]
        else:
            alignment = msa['alignment']

        # create new array for confidence
        confidence_matrix = []
        character_matrix = []

        # iterate over each taxon
        for i, taxon in enumerate(taxa):
            idx = alms.taxa.index(taxon) + 1

            # get the numerical sequence
            nums = alignment[i]

            # store confidences per line
            confidences = []

            # store chars per line
            chars = []

            # iterate over the sequence
            for j, num in enumerate(nums):
                col = [alm[j] for alm in alignment]
                score = 0
                count = 0

                # get the char
                if num != '-':
                    charA = dotjoin(taxa[i], msa['alignment'][i][j], num.split('.')[2])
                    chars += [charA]
                    try:
                        occs[charA] += [concept]
                    except:
                        occs[charA] = [concept]
                else:
                    chars += ['-']

                for k, numB in enumerate(col):
                    if k != i:
                        if num == '-' and numB == '-':
                            pass
                        else:
                            if numB != '-' and num != '-':
                                # get the second char
                                charB = dotjoin(
                                    taxa[k], msa['alignment'][k][j], numB.split('.')[2])
                                try:
                                    corrs[charA][charB] += 1
                                except:
                                    try:
                                        corrs[charA][charB] = 1
                                    except:
                                        corrs[charA] = {charB: 1}

                            gaps = False
                            if num == '-' and numB != '-':
                                numA = charstring(idx)
                                gaps = True
                            elif numB == '-' and num != '-':
                                numB = charstring(alms.taxa.index(taxa[k]))
                                numA = num
                                gaps = True
                            else:
                                numA = num

                            scoreA = scorer[numA, numB]
                            scoreB = scorer[numB, numA]
                            this_score = max(scoreA, scoreB)

                            if not gaps:
                                score += this_score
                                count += 1
                            else:
                                score += this_score * gap_weight
                                count += gap_weight

                if count:
                    score = score / count
                else:
                    score = -25

                confidences += [int(score + 0.5)]
                values += [int(score + 0.5)]
            confidence_matrix += [confidences]
            character_matrix += [chars]

        # append confidence matrix to alignments
        alms.msa[ref][key]['confidence'] = confidence_matrix
        alms.msa[ref][key]['_charmat'] = character_matrix

    # sort the values
    values = sorted(set(values + [1]))

    # make conversion to scale of 100 values
    converter = {}
    valsA = values[:values.index(1)]
    valsB = values[values.index(1):]
    stepA = 50 / (len(valsA) + 1)
    stepB = 75 / (len(valsB) + 1)
    for i, score in enumerate(valsA):  # values[:values.index(0)):
        converter[score] = int((stepA * i) / 4 + 0.5)
    for i, score in enumerate(valsB):
        converter[score] = int(stepB * i + 0.5) + 50

    # iterate over keys again
    for key, msa in alms.msa[ref].items():
        # get basic stuff
        for i, line in enumerate(msa['confidence']):
            for j, cell in enumerate(line):
                alms.msa[ref][key]['confidence'][i][j] = converter[cell]

    jsond = {}
    for key, corr in corrs.items():
        splits = [c.split('.') + [o] for c, o in corr.items()]
        sorts = sorted(splits, key=lambda x: (x[0], -x[3]))
        new_sorts = []

        # check for rowspan
        spans = {}
        for a, b, c, d in sorts:
            if a in spans:
                if spans[a] < 3 and d > 1:
                    spans[a] += 1
                    new_sorts += [[a, b, c, d]]
            else:
                if d > 1:
                    spans[a] = 1
                    new_sorts += [[a, b, c, d]]

        bestis = []
        old_lang = ''
        counter = 0
        for a, b, c, d in new_sorts:
            new_lang = a
            if new_lang != old_lang:
                old_lang = new_lang

                tmp = '<tr class="display">'
                tmp += '<td class="display" rowspan={0}>'.format(spans[a])
                tmp += a + '</td>'
                tmp += '<td class="display" onclick="show({0});"><span '.format(
                    "'" + dotjoin(a, b, c) + "'")
                tmp += 'class="char {0}">' + b + '</span></td>'
                tmp += '<td class="display">'
                tmp += c + '</td>'
                tmp += '<td class="display">' + str(d) + '</td>'
                tmp += '<td class="display">' + str(len(occs[dotjoin(a, b, c)])) + '</td>'
                tmp += '</tr>'
                t = 'dolgo_' + token2class(b, rcParams['dolgo'])

                # bad check for three classes named differently
                if t == 'dolgo__':
                    t = 'dolgo_X'
                elif t == 'dolgo_1':
                    t = 'dolgo_TONE'
                elif t == 'dolgo_0':
                    t = 'dolgo_ERROR'

                bestis += [tmp.format(t)]
                counter += 1

            elif counter > 0:
                tmp = '<tr class="display">'
                tmp += '<td class="display" onclick="show({0});"><span '.format(
                    "'" + dotjoin(a, b, c) + "'")
                tmp += 'class="char {0}">' + b + '</span></td>'
                tmp += '<td class="display">' + c + '</td>'
                tmp += '<td class="display">' + str(d) + '</td>'
                tmp += '<td class="display">' + str(len(occs[dotjoin(a, b, c)])) + '</td>'
                tmp += '</tr>'

                t = 'dolgo_' + token2class(b, rcParams['dolgo'])

                # bad check for three classes named differently
                if t == 'dolgo__':
                    t = 'dolgo_X'
                elif t == 'dolgo_1':
                    t = 'dolgo_TONE'
                elif t == 'dolgo_0':
                    t = 'dolgo_ERROR'

                bestis += [tmp.format(t)]
                counter += 1
                old_lang = new_lang
            else:
                old_lang = new_lang
                counter = 0

        jsond[key] = [''.join(bestis), occs[key]]

    return jsond


def get_correspondences(alms, ref='lexstatid'):
    """
    Compute sound correspondences for a given set of aligned cognates.
    """
    # store all correspondences
    corrs = {}

    # store occurrences
    occs = {}

    for key, msa in alms.msa[ref].items():
        # get basic stuff
        idxs = msa['ID']
        taxa = msa['taxa']
        concept = cgi.escape(alms[idxs[0], 'concept'], True)

        # get numerical representation of alignments
        if 'numbers' in alms.header:
            alignment = [class2tokens(
                alms[idxs[i], 'numbers'],
                msa['alignment'][i]) for i in range(len(idxs))]
        else:
            alignment = msa['alignment']

        # create new array for confidence
        character_matrix = []

        # iterate over each taxon
        for i, taxon in enumerate(taxa):
            # get the numerical sequence
            nums = alignment[i]

            # store chars per line
            chars = []

            # iterate over the sequence
            for j, num in enumerate(nums):
                col = [alm[j] for alm in alignment]

                # get the char
                if num != '-':
                    charA = dotjoin(taxa[i], msa['alignment'][i][j], num.split('.')[2])
                    chars += [charA]
                    try:
                        occs[charA] += [concept]
                    except:
                        occs[charA] = [concept]
                else:
                    chars += ['-']

                for k, numB in enumerate(col):
                    if k != i:
                        if num == '-' and numB == '-':
                            pass
                        else:
                            if numB != '-' and num != '-':
                                # get the second char
                                charB = dotjoin(
                                    taxa[k],
                                    msa['alignment'][k][j],
                                    numB.split('.')[2])
                                try:
                                    corrs[charA][charB] += 1
                                except:
                                    try:
                                        corrs[charA][charB] = 1
                                    except:
                                        corrs[charA] = {charB: 1}

            character_matrix += [chars]

        # append confidence matrix to alignments
        alms.msa[ref][key]['_charmat'] = character_matrix

    return corrs, occs
