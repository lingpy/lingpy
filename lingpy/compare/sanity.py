"""
Module provides basic checks for wordlists.
"""
from __future__ import (
        unicode_literals, print_function, absolute_import, division)
import networkx as nx
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.approximation.clique import max_clique
from itertools import combinations
from collections import defaultdict
from lingpy import log


def _mutual_coverage(taxA, taxB, wordlist, concepts):
    return set(
        [w for w in wordlist.get_list(
            col=taxA, 
            flat=True, 
            entry=concepts
            ) if w in wordlist.get_list(
                col=taxB, 
                flat=True,
                entry=concepts)])

def _get_concepts(wordlist, concepts):
    return {c: set(wordlist.get_list(col=c, flat=True, entry=concepts)) for c in
        wordlist.cols}

def mutual_coverage(wordlist, concepts='concept'):
    """Compute mutual coverage for all language pairs in your data.
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        Your Wordlist object (or a descendant class).
    concepts : str (default="concept")
        The column which stores your concepts.

    Returns
    -------
    coverage : dict
        A dictionary of dictionaries whose value is the number of items two
        languages share.

    Examples
    --------
    
    Compute coverage for the KSL.qlc dataset::
      
      >>> from lingpy.compare.sanity import mutual_coverage
      >>> from lingpy import *
      >>> from lingpy.tests.util import test_data
      >>> wl = Wordlist(test_data('KSL.qlc'))
      >>> cov = mutual_coverage(wl)
      >>> cov['English']['German']
      200

    See also
    --------
    mutual_coverage_check
    mutual_coverage_subset
    average_coverage
    """
    coverage = defaultdict(dict)
    concepts = _get_concepts(wordlist, concepts)
    for t1, t2 in combinations(wordlist.cols, r=2):
        coverage[t1][t2] = concepts[t1].intersection(concepts[t2])
        coverage[t2][t1] = coverage[t1][t2]
    return coverage

def mutual_coverage_check(wordlist, threshold, concepts='concept'):
    """Check whether a given mutual coverage is fulfilled by the dataset.
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        Your Wordlist object (or a descendant class).
    concepts : str (default="concept")
        The column which stores your concepts.
    threshold : int
        The threshold which should be checked.
    
    Returns
    -------
        c: bool
            True, if coverage is fulfilled for all language pairs, False if
            otherwise.

    Examples
    --------
    Compute minimal mutual coverage for the KSL dataset::

      >>> from lingpy.compare.sanity import mutual_coverage
      >>> from lingpy import *
      >>> from lingpy.tests.util import test_data
      >>> wl = Wordlist(test_data('KSL.qlc'))
      >>> for i in range(wl.height, 1, -1):
              if mutual_coverage_check(wl, i):
                  print('mutual coverage is {0}'.format(i))
                  break
          200

    See also
    --------
    mutual_coverage
    mutual_coverage_subset
    average_coverage
    """
    mc = mutual_coverage(wordlist, concepts)
    for coverage in mc.values():
        if [x for x in coverage if len(coverage[x]) < threshold]:
            return False
    return True

def mutual_coverage_subset(wordlist, threshold, concepts='concept'):
    """Compute maximal mutual coverage for all language in a wordlist.
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        Your Wordlist object (or a descendant class).
    concepts : str (default="concept")
        The column which stores your concepts.
    threshold : int
        The threshold which should be checked.  

    Returns
    -------
    coverage : tuple
        A tuple consisting of the number of languages for which the coverage
        could be found as well as a list of all pairings in which this coverage
        is possible. The list itself contains the mutual coverage inside each
        pair and the list of languages.

    Examples
    --------
    Compute all sets of languages with coverage at 200 for the KSL dataset::

      >>> from lingpy.compare.sanity import mutual_coverage_subset
      >>> from lingpy import *
      >>> from lingpy.tests.util import test_data
      >>> wl = Wordlist(test_data('KSL.qlc'))
      >>> number_of_languages, pairs = mutual_coverage_subset(wl, 200)
      >>> for number_of_items, languages in pairs:
              print(number_of_items, ','.join(languages))
          200 Albanian,English,French,German,Hawaiian,Navajo,Turkish

    See also
    --------
    mutual_coverage
    mutual_coverage_check
    average_coverage
    """
    coverage = mutual_coverage(wordlist, concepts)

    G = nx.Graph()
    for tax in wordlist.cols:
        G.add_node(tax)
    for taxA, taxB in combinations(wordlist.cols, r=2):
        if len(coverage[taxA][taxB]) >= threshold:
            G.add_edge(taxA, taxB, coverage=coverage[taxA][taxB])
    
    best_cliques = defaultdict(list)
    best_clique = 0
    for clique in find_cliques(G):
        sums = []
        for taxA, taxB in combinations(clique, r=2):
            sums += [len(G[taxA][taxB]['coverage'])]
        if sums:
            val = int(sum(sums) / len(sums) + 0.5)
            best_cliques[len(clique)] += [(val, sorted(clique))]
            if len(clique) > best_clique:
                best_clique = len(clique)
    return best_clique, best_cliques[best_clique]


def average_coverage(wordlist, concepts='concepts'):
    """Compute average mutual coverage for a given wordlist.
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        Your Wordlist object (or a descendant class).
    concepts : str (default="concept")
        The column which stores your concepts.

    Returns
    -------
    coverage : dict
        A dictionary of dictionaries whose value is the number of items two
        languages share.

    Examples
    --------
    
    Compute coverage for the KSL.qlc dataset::
      
      >>> from lingpy.compare.sanity import average_coverage
      >>> from lingpy import *
      >>> from lingpy.tests.util import test_data
      >>> wl = Wordlist(test_data('KSL.qlc'))
      >>> average_coverage(wl)
      1.0

    See also
    --------
    mutual_coverage_check
    mutual_coverage_subset
    mutual_coverage

    """
    mc = mutual_coverage(wordlist)
    score = []
    for v in mc.values():
        for key, val in v.items():
            score += [len(val)]
    return sum(score) / len(score) / wordlist.height


def synonymy(wordlist, concepts='concept', languages='doculect'):
    """Check the number of synonyms per language and concept.
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        Your Wordlist object (or a descendant class).
    concepts : str (default="concept")
        The column which stores your concepts.
    languages : str (default="doculect")
        The column which stores your language names. 
    
    Returns
    -------
    synonyms : dict
        A dictionary with language and concept as key and the number of
        synonyms as value.

    Examples
    --------
    Calculate synonymy in KSL.qlc dataset::

      >>> from lingpy.compare.sanity import synonymy      
      >>> from lingpy import *
      >>> from lingpy.tests.util import test_data
      >>> wl = Wordlist(test_data('KSL.qlc'))
      >>> syns = synonymy(wl)
      >>> for a, b in syns.items():
              if b > 1:
                  print(a[0], a[1], b)

    There is no case where synonymy exceeds 1 word per concept per language,
    since :evobib:`Kessler2001` was paying particular attention to avoid
    synonyms.
    """
    synonyms = defaultdict(int)
    for idx, language, concept in wordlist.iter_rows(languages, concepts):
        synonyms[language, concept] += 1

    return synonyms


def check_length(a, b, dimA=1, dimB=1):
    """
    Custom function to check the length of two basictypes in LingPy.
    """
    if dimA < 3 or dimB < 3:
        entityA, entityB = a if dimA == 1 else a.n, b if dimB == 1 else b.n
        if len(entityA) != len(entityB):
            return False
    else:
        if len(a) == len(b):
            if len(a.n) == len(b.n):
                for pA, pB in zip(a.n, b.n):
                    if len(pA) != len(pB):
                        return False
            else:
                return False
        else:
            return False

    return True
                        


def check_sequence_length(
        wordlist, 
        entities=['tokens', 'crossids', 'morphemes', 'structure'],
        dimensions=[2, 1, 2, 1]
        ):
    """Function checks for identical sequence length in different columns.
    """
    fails, errors = [], 0
    for (eA, dA), (eB, dB) in combinations(zip(entities, dimensions), r=2):
        for idx in wordlist:
            if not check_length(
                    wordlist[idx, eA], 
                    wordlist[idx, eB],
                    dA,
                    dB
                    ):
                errors += 1
                log.warning(
                        '{0} | {1} | {2} | {3} | {4} | {5}'.format(
                            errors,
                            idx,
                            eA,
                            eB,
                            wordlist[idx, eA],
                            wordlist[idx, eB]
                            )
                        )
                fails += [idx]
    return fails


def check_cognates(wordlist, ref='crossids'):
    """Function checks for internal consistency of partial cognates."""
    fails = []
    for idx, cogids in wordlist.iter_rows(ref):
        if len(set(cogids)) != len(cogids):
            log.warning('duplicates in {0}'.format(cogids))
            fails += [idx]
    return fails


def check_strict_cognates(
        wordlist,
        ref='crossids',
        segments='tokens'
        ):
    """Check if cognates are really strict."""
    fails, errors = [], 0
    etd = wordlist.get_etymdict(ref=ref)
    for cogid in etd:
        for vals in etd[cogid]:
            if vals:
                if not str(cogid).isdigit():
                    fails += vals
                else:
                    alms = []
                    for idx in vals:
                        try:
                            tokens = wordlist[idx, segments].n[wordlist[idx,
                                ref].index(cogid)]
                            alms += [(idx, tokens)]
                        except:
                            fails += [idx]
                    for idx, tokens in alms[1:]:
                        if str(tokens) != str(alms[0][1]):
                            fails += [idx]
                            errors += 1
                            log.warning('{0} | {1} | {2} | {3:15} | {4:15}'.format(
                                        errors, idx, alms[0][0],
                                        str(tokens), str(alms[0][1])))
    return fails

