"""
Module provides methods for the evaluation of automatic linguistic reconstruction analyses.
"""
from __future__ import unicode_literals, division
import logging

from ..settings import rcParams
from ..align.pairwise import edit_dist
from ..sequence.sound_classes import ipa2tokens, tokens2class
from .. import log
from ..util import setdefaults, as_string


def mean_edit_distance(
        wordlist,
        gold="proto",
        test="consensus",
        ref="cogid",
        tokens=True,
        classes=False,
        **keywords):
    """
    Function computes the edit distance between gold standard and test set.

    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        The wordlist object containing the data for a given analysis.
    gold : str (default="proto")
        The name of the column containing the gold-standard solutions.
    test = "consensus"
        The name of the column containing the test solutions.
    
    Returns
    -------
    dist : float
        The mean edit distance between gold and test reconstructions.

    Notes
    -----
    This function has an alias ("med"). Calling it will produce the same
    results.
    """
    setdefaults(
        keywords,
        merge_vowels=rcParams['merge_vowels'],
        model=rcParams['model'])

    distances = []

    for key, idxs in wordlist.get_etymdict(ref=ref).items():
        # get only valid numbers for index-search
        idx = [idx[0] for idx in idxs if idx != 0][0]

        log.debug('{0}, {1}'.format(idx, idxs))

        # get proto and consensus from wordlist
        proto = wordlist[idx, gold]
        consensus = wordlist[idx, test]

        log.debug('{0}, {1}'.format(proto, consensus))

        if tokens or classes:
            proto = ipa2tokens(proto, **keywords)
            consensus = ipa2tokens(consensus, **keywords)

            if classes:
                proto = tokens2class(proto, **keywords)
                consensus = tokens2class(consensus, **keywords)

        distances.append(edit_dist(proto, consensus, normalized=False))

    med = sum(distances) / len(distances)
    log.info('MEAN ED: {0:.2f}'.format(med))
    return med


def med(wordlist, **keywords):
    # alias for mean_edit_distance
    return mean_edit_distance(wordlist, **keywords)
