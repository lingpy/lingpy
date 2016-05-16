# *-* coding: utf-8 *-*
"""
Evaluation methods for automatic cognate detection.
"""
from __future__ import print_function, division, unicode_literals
import codecs
from itertools import combinations
from collections import defaultdict

import logging
from lingpy import log
from lingpy.util import identity, as_string 

def _get_bcubed_score(one, other):
    tmp = defaultdict(list)
    for x, y in zip(one, other):
        tmp[x].append(y)
    bcp = 0.0
    for x in tmp:
        for y in tmp[x]:
            bcp += tmp[x].count(y) / len(tmp[x])
    return bcp / len(other)

def _get_cogs(ref, concept, modify_ref, wordlist):
    idxs = wordlist.get_list(row=concept, flat=True)
    bidx = [i + 1 for i in range(len(idxs))]
    cogs = wordlist.get_list(row=concept, entry=ref, flat=True)
    tmp = {}
    for a, b in zip(cogs, bidx):
        if modify_ref(a) not in tmp:
            tmp[(a)] = b
    return [tmp[modify_ref(i)] for i in cogs]

def _format_results(results, p, r, f):
    """
    Print out the results of an analysis.
    """

    return """*************************')
* {0:7}-Scores        *
* --------------------- *
* Precision:     {1:.4f} *
* Recall:        {2:.4f} *
* F-Scores:      {3:.4f} *
*************************'""".format(
        results, p, r, f)

def bcubes(wordlist, gold='cogid', test='lexstatid', modify_ref=False, pprint=True, 
        per_concept=False):
    """
    Compute B-Cubed scores for test and reference datasets.

    Parameters
    ----------
    lex : :py:class:`lingpy.basic.wordlist.Wordlist`
        A :py:class:`lingpy.basic.wordlist.Wordlist` class or a daughter class,
        (like the :py:class:`~lingpy.compare.lexstat.LexStat` class used for the
        computation). It should have two columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    modify_ref : function (default=False)
        Use a function to modify the reference. If your cognate identifiers
        are numerical, for example, and negative values are assigned as
        loans, but you want to suppress this behaviour, just set this
        keyword to "abs", and all cognate IDs will be converted to their
        absolute value.
    pprint : bool (default=True)
        Print out the results
    per_concept : bool (default=False)
        Compute b-cubed scores per concep and not for the whole data in one
        piece.

    Returns
    -------
    t : tuple
        A tuple consisting of the precision, the recall, and the harmonic mean
        (F-scores).

    Notes
    -----
    B-Cubed scores were first described by :evobib:`Bagga1998` as part of an
    algorithm. Later on, :evobib:`Amigo2009` showed that they can also used as
    to compare cluster decisions. :evobib:`Hauer2011` applied the B-Cubed
    scores first to the task of automatic cognate detection.
    
    See also
    --------
    diff
    pairs
    """
    # if loans are treated as homologs
    evl = modify_ref if modify_ref else identity
    
    def get_scores(one, other):
        for _, line in wordlist.get_etymdict(ref=one, modify_ref=modify_ref).items():
            line = [value for value in [evl(x[0]) for x in line if x != 0]]
            # check for linesize
            if len(line) > 1:
                # get cognate-ids in the other set for the line
                other_line = [evl(wordlist[idx, other]) for idx in line]

                # get the recall
                for idx in other_line:
                    yield other_line.count(idx) / len(line)
            else:
                yield 1.0

    if per_concept:
        bcr, bcp, fsc = [], [], []
        for concept in wordlist.rows:
            idxsG = _get_cogs(gold, concept, evl, wordlist)
            idxsT = _get_cogs(test, concept, evl, wordlist)
            r = _get_bcubed_score(idxsG, idxsT)
            p = _get_bcubed_score(idxsT, idxsG)
            f = 2 * ((r * p) / (p + r))
            bcr += [r]
            bcp += [p]
            fsc += [f]
            
            as_string('{0:15}\t{1:.2f}\t{2:.2f}\t{3:.2f}'.format(
                    concept, p, r, f), pprint=pprint)
    else:
        # b-cubed recall
        bcr = list(get_scores(gold, test))
        # b-cubed precision
        bcp = list(get_scores(test, gold))
        fsc = []

    # calculate general scores
    BCP = sum(bcp) / len(bcp)
    BCR = sum(bcr) / len(bcr)
    FSC = sum(fsc) / len(fsc) if fsc else 2 * ((BCP * BCR) / (BCP + BCR))
    
    as_string(_format_results('B-Cubed', BCP, BCR, FSC), pprint=pprint)

    return BCP, BCR, FSC

def partial_bcubes(wordlist, gold, test, pprint=True):
    """
    Compute B-Cubed scores for test and reference datasets for partial cognate\
            detection.

    Parameters
    ----------
    wordlist : :py:class:`~lingpy.basic.wordlist.Wordlist`
        A :py:class:`~lingpy.basic.wordlist.Wordlist`, or one of it's daughter
        classes (like, e.g., the :py:class:`~lingpy.compare.partial.Partial`
        class used for computation of partial cognates. It should have two
        columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    pprint : bool (default=True)
        Print out the results

    Returns
    -------
    t : tuple
        A tuple consisting of the precision, the recall, and the harmonic mean
        (F-scores).

    Notes
    -----
    B-Cubed scores were first described by :evobib:`Bagga1998` as part of an
    algorithm. Later on, :evobib:`Amigo2009` showed that they can also used as
    to compare cluster decisions. :evobib:`Hauer2011` applied the B-Cubed
    scores first to the task of automatic cognate detection.
    
    See also
    --------
    bcubes
    diff
    pairs
    """
    
    # here's the point with bcubes for fuzzy: if we compare, we need to make
    # sure we count whether one instance is identical, not whether all of them
    # are identical!
    
    def get_scores(one, other):
        scores = []
        multiple_items = []
        for k,v in wordlist.get_etymdict(ref=one).items():
            _idxs = [val for val in v if val != 0]
            # now we need to get the position in the index
            poss,idxs = [],[]
            for val in _idxs:
                if len(val) > 1:
                    multiple_items += [len(val)]
                for idx in val:
                    new_pos = [i for i,cog in zip(range(len(wordlist[idx,one])),
                        wordlist[idx,one]) if cog == k]
                    idxs += [idx for x in new_pos]
                    poss += new_pos
            if len(idxs) > 1:
                other_idxs = [wordlist[idx,other][pos] for pos,idx in zip(poss,idxs)]
                for idx in other_idxs:
                    scores += [other_idxs.count(idx) / len(idxs)]
            else: 
                scores += [1]
        return sum(scores) / len(scores)

    bcr = get_scores(gold, test)
    bcp = get_scores(test, gold)
    bcf = 2 * ((bcp * bcr) / (bcp + bcr))
    
    as_string(_format_results('B-Cubed', bcp, bcr, bcf), 
            pprint=pprint)
    return bcp, bcr, bcf

def pairs(lex, gold='cogid', test='lexstatid', modify_ref=False, pprint=True,
        _return_string=False):
    """
    Compute pair scores for the evaluation of cognate detection algorithms.
    
    Parameters
    ----------
    lex : :py:class:`lingpy.compare.lexstat.LexStat`
        The :py:class:`~lingpy.compare.lexstat.LexStat` class used for the
        computation. It should have two columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    modify_ref : function (default=False)
        Use a function to modify the reference. If your cognate identifiers
        are numerical, for example, and negative values are assigned as
        loans, but you want to suppress this behaviour, just set this
        keyword to "abs", and all cognate IDs will be converted to their
        absolute value.
    pprint : bool (default=True)
        Print out the results

    Returns
    -------
    t : tuple
        A tuple consisting of the precision, the recall, and the harmonic mean
        (F-scores).
    
    Notes
    -----
    Pair-scores can be computed in different ways, with often different
    results. This variant follows the description by :evobib:`Bouchard-Cote2013`.
    
    See also
    --------
    diff
    bcubes
    """
    # if loans are treated as homologs
    evl = modify_ref if modify_ref else identity

    def get_pairs(ref):
        for key, line in lex.get_etymdict(ref=ref, modify_ref=modify_ref).items():
            line = [value for value in [evl(x[0]) for x in line if x != 0]]
            for a, b in combinations(line, r=2):
                yield tuple(sorted([a, b]))

    pairsG = set(get_pairs(gold))
    pairsT = set(get_pairs(test))
    
    # calculate precision and recall
    pp = len(pairsG.intersection(pairsT)) / len(pairsT)
    pr = len(pairsG.intersection(pairsT)) / len(pairsG)
    fs = 2 * (pp * pr) / (pp + pr)

    # print the results if this option is chosen
    as_string(_format_results('Pairs', pp, pr, fs), pprint=pprint)
    
    return pp, pr, fs


def diff(
        wordlist,
        gold='cogid',
        test='lexstatid',
        modify_ref=False,
        pprint=True,
        filename='',
        tofile=True,
        transcription="ipa"):
    r"""
    Write differences in classifications on an item-basis to file.

    lex : :py:class:`lingpy.compare.lexstat.LexStat`
        The :py:class:`~lingpy.compare.lexstat.LexStat` class used for the
        computation. It should have two columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    modify_ref : function (default=False)
        Use a function to modify the reference. If your cognate identifiers
        are numerical, for example, and negative values are assigned as
        loans, but you want to suppress this behaviour, just set this
        keyword to "abs", and all cognate IDs will be converted to their
        absolute value.
    pprint : bool (default=True)
        Print out the results
    filename : str (default='')
        Name of the output file. If not specified, it is identical with the
        name of the :py:class:`~lingpy.compare.lexstat.LexStat`, but with the
        extension ``diff``.
    tofile : bool (default=True)
        If set to c{False}, no data will be written to file, but instead, the
        data will be returned.
    transcription : str (default="ipa")
        The file in which the transcriptions are located (should be a string,
        no segmentized version, for convenience of writing to file).

    Returns
    -------
    t : tuple
        A nested tuple consisting of two further tuples. The first
        containing precision, recall, and harmonic mean
        (F-scores), the second containing the same values for the pair-scores.

    Notes
    -----
    If the **tofile** option is chosen, the results are written to a specific
    file with the extension ``diff``. This file contains all cognate sets in
    which there are differences between gold standard and test sets. It also
    gives detailed information regarding false positives, false negatives, and
    the words involved in these wrong decisions.

    .. This function also calculates the "transformation" score. This score is
    .. based on the calculation of steps that are needed to transform one cluster
    .. for one set of meanings into the other. Ideally, if there are *n* different
    .. cognate sets covering one gloss in the gold standard, the minimal length of
    .. a mapping to convert the *m* cognate sets of the test set into the gold standard
    .. is *n*. In this case, both gold standard and test set are identical.
    .. However, if gold standard and test set differ, the number of mappings
    .. necessarily exceeds *m* and *n*. Based on this, the transformation
    .. precision is defined as :math:`\frac{m}{M}`, where *m* is the number of
    .. distinct clusters in the test set and *M* is the length of the mapping.
    .. Accordingly, the recall is defined as :math:`\frac{n}{M}`, where *n* is the
    .. number of clusters in the gold standard.

    .. Note that if precision is lower than 1.0, this means there are false
    .. positive decisions in the test set. Accordingly, a recall lower than 1.0
    .. indicates that there are false negative decisions in the test set.
    .. The drawback of this score is that it is not sensitive regarding the
    .. distinct number of decisions in which gold standard and test set differ, so
    .. the recall can be very low although most of the words have been grouped
    .. accurately. The advantage is that it can be directly interpreted in terms
    .. of 'false positive/false negative' decisions.

    See also
    --------
    bcubes
    pairs
    """
    filename = filename or wordlist.filename
    loan = modify_ref if modify_ref else identity

    # open file
    if tofile:
        f = codecs.open(filename + '.diff', 'w', 'utf-8')

    # get a formatter for language names
    lform = '{0:' + str(max([len(l) for l in wordlist.cols])) + '}'
    
    preT, recT = [], []
    preB, recB = [], []
    preP, recP = [], []

    def get_pairs(cogs, idxs):
        tmp = defaultdict(list)
        for x, y in zip(cogs, idxs):
            tmp[x].append(y)
        for x in tmp:
            for yA, yB in combinations(tmp[x], r=2):
                yield tuple(sorted([yA, yB]))

    for concept in wordlist.rows:
        idxs = wordlist.get_list(row=concept, flat=True)
        # get the basic index for all seqs
        bidx = [i + 1 for i in range(len(idxs))]

        cogsG = _get_cogs(gold, concept, loan, wordlist)
        cogsT = _get_cogs(test, concept, loan, wordlist)

        if cogsG != cogsT:
            # calculate the transformation distance of the sets
            tramGT = len(set(zip(cogsG, cogsT)))
            tramG = len(set(cogsG))
            tramT = len(set(cogsT))
            preT += [tramT / tramGT]
            recT += [tramG / tramGT]

            # calculate the bcubed precision for the sets
            preB += [_get_bcubed_score(cogsT, cogsG)]

            # calculate b-cubed recall
            recB += [_get_bcubed_score(cogsG, cogsT)]

            # calculate pair precision
            pairsG = set(get_pairs(cogsG, idxs))
            pairsT = set(get_pairs(cogsT, idxs))

            preP.append(len(pairsT.intersection(pairsG)) / len(pairsT) if pairsT else 1.0)
            recP.append(len(pairsT.intersection(pairsG)) / len(pairsG) if pairsG else 1.0)
            fp = "no" if preP[-1] == 1.0 else "yes"
            fn = "no" if recP[-1] == 1.0 else "yes"

            if tofile:
                f.write(
                    "Concept: {0}, False Positives: {1}, False Negatives: {2}\n".format(
                        concept, fp, fn))

            # get the words
            words = [wordlist[i, 'ipa'] for i in idxs]
            langs = [wordlist[i, 'taxa'] for i in idxs]

            # get a word-formater
            wform = '{0:' + str(max([len(w) for w in words])) + '}'

            # write differences to file
            if tofile:
                for word, lang, cG, cT in sorted(
                        zip(words, langs, cogsG, cogsT),
                        key=lambda x: (x[2], x[3])):
                    f.write('{0}\t{1}\t{2:4}\t{3:4}\n'.format(
                        lform.format(lang), wform.format(word), cG, cT))
                f.write('#\n')
        else:
            preT += [1.0]
            recT += [1.0]
            preB += [1.0]
            recB += [1.0]
            preP += [1.0]
            recP += [1.0]

    bp = sum(preB) / len(preB)
    br = sum(recB) / len(recB)
    bf = 2 * (bp * br) / (bp + br)
    pp = sum(preP) / len(preP)
    pr = sum(recP) / len(recP)
    pf = 2 * (pp * pr) / (pp + pr)
    

    as_string(_format_results('B-Cubed', bp, br, bf) + \
            _format_results('Pair', pp, pr, pf), 
            pprint=pprint)

    if tofile:
        f.write('B-Cubed Scores:\n')
        f.write('Precision: {0:.4f}\n'.format(bp))
        f.write('Recall:    {0:.4f}\n'.format(br))
        f.write('F-Score:   {0:.4f}\n'.format(bf))
        f.write('#\n')
        f.write('Pair Scores:\n')
        f.write('Precision: {0:.4f}\n'.format(pp))
        f.write('Recall:    {0:.4f}\n'.format(pr))
        f.write('F-Score:   {0:.4f}\n'.format(pf))
        f.close()
        log.file_written(filename + '.diff')
    else:
        return (bp, br, bf), (pp, pr, pf)

def npoint_ap(scores, cognates, reverse=False):
    """
    Calculate the n-point average precision.
    
    Parameters
    ----------
    scores : list
        The scores of your algorithm for pairwise string comparison. 
    cognates : list
        The cognate codings of the word pairs you compared. 1 indicates that
        the pair is cognate, 0 indicates that it is not cognate.
    reverse : bool (default=False)
        The order of your ranking mechanism. If your algorithm yields high
        scores for words which are probably cognate, and low scores for
        non-cognate words, you should set this keyword to "True".

    Notes
    -----
    This follows the description in :evobib:`Kondrak2002`. The n-point average
    precision is useful to compare the discriminative force of different
    algorithms for string similarity, or to train the parameters of a given
    algorithm.

    Examples
    --------
    
    >>> scores = [1, 2, 3, 4, 5]
    >>> cognates = [1, 1, 1, 0, 0]
    >>> from lingpy.evaluate.acd import npoint_ap
    >>> npoint_ap(scores, cognates)
    1.0

    """
    p = 0.0
    cognate_count = 0
    for k,(score, cognate) in enumerate(sorted(zip(scores, cognates),
        key=lambda x: x[0], reverse=reverse)):
        if cognate == 1:
            cognate_count += 1
            p += cognate_count / (k+1.0)
    try: 
        return p / cognates.count(1)
    except ZeroDivisionError:
        log.warn("Encountered Zero Division in npoint_ap, your data seems to contain no cognates.")
        return 0

