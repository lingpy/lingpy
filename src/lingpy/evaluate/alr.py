# *-* coding: utf-8 *-*
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
        
    stress : str (default=rcParams['stress'])
        A string containing the stress symbols used in the sound-class
        conversion. Defaults to the stress as defined in
        ~lingpy.settings.rcParams.

    diacritics : str (default=rcParams['diacritics'])
        A string containing diacritic symbols used in the sound-class
        conversion. Defaults to the diacritic symbolds defined in
        ~lingpy.settings.rcParams.

    cldf : bool (default=False)
        If set to True, this will allow for a specific treatment of phonetic
        symbols which cannot be completely resolved (e.g., laryngeal h₂ in
        Indo-European). Following the `CLDF <http://cldf.clld.org>`_
        specifications (in particular the specifications for writing
        transcriptions in segmented strings, as employed by the `CLTS
        <http://calc.digling.org/clts/>`_ initiative), in cases of insecurity
        of pronunciation, users can adopt a ```source/target``` style, where
        the source is the symbol used, e.g., in a reconstruction system, and
        the target is a proposed phonetic interpretation. This practice is also
        accepted by the `EDICTOR <http://edictor.digling.org>`_ tool.

    
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
        model=rcParams['model'],
        stress=rcParams['stress'],
        diacritics=rcParams['diacritics'],
        cldf=False)

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
                proto = tokens2class(proto, keywords['model'], stress=keywords['stress'],
                        diacritics=keywords['diacritics'],
                        cldf=keywords['cldf'])
                consensus = tokens2class(consensus, keywords['model'],
                        stress=keywords['stress'],
                        diacritics=keywords['diacritics'],
                        cldf=keywords['cldf'])

        distances.append(edit_dist(proto, consensus, normalized=False))

    med = sum(distances) / len(distances)
    log.info('MEAN ED: {0:.2f}'.format(med))
    return med


def med(wordlist, **keywords):
    # alias for mean_edit_distance
    return mean_edit_distance(wordlist, **keywords)
