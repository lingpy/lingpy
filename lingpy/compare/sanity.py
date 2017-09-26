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
    coverage = defaultdict(dict)
    concepts = _get_concepts(wordlist, concepts)
    for t1, t2 in combinations(wordlist.cols, r=2):
        coverage[t1][t2] = len(concepts[t1].intersection(concepts[t2]))
        coverage[t2][t1] = coverage[t1][t2]
    return coverage

def mutual_coverage_check(wordlist, threshold, concepts='concept'):
    """Check whether a given mutual coverage is fulfilled by the dataset.

    Returns
    -------
        c: bool
            True, if coverage is fulfilled for all language pairs, False if
            otherwise.
    """
    mc = mutual_coverage(wordlist, concepts)
    for coverage in mc.values():
        if [x for x in coverage if coverage[x] < threshold]:
            return False
    return True

def mutual_coverage_subset(wordlist, threshold, concepts='concept'):
    """Compute maximal mutual coverage for all language in a wordlist.
    
    Note
    ----
    Returns all languages in a sample for which coverage is minimally as
    defined in the threshold. Coverage means the number of concepts for which
    there is a translation in both language A and language B.
    """
    coverage = mutual_coverage(wordlist, concepts)

    G = nx.Graph()
    for tax in wordlist.cols:
        G.add_node(tax)
    for taxA, taxB in combinations(wordlist.cols, r=2):
        if coverage[taxA][taxB] >= threshold:
            G.add_edge(taxA, taxB, coverage=coverage[taxA][taxB])
    
    best_cliques = defaultdict(list)
    best_clique = 0
    for clique in find_cliques(G):
        sums = []
        for taxA, taxB in combinations(clique, r=2):
            sums += [G.edge[taxA][taxB]['coverage']]
        if sums:
            val = int(sum(sums) / len(sums) + 0.5)
            best_cliques[len(clique)] += [(val, sorted(clique))]
            if len(clique) > best_clique:
                best_clique = len(clique)
    return best_clique, best_cliques[best_clique]

def synonymy(wordlist, segments='tokens', sound_classes=False):
    words = defaultdict(list)
    pass
