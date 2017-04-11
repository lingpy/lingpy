"""
Module provides basic checks for wordlists.
"""
from __future__ import (
        unicode_literals, print_function, absolute_import, division)
import networkx as nx
from networkx.algorithms.clique import find_cliques
from itertools import combinations
from collections import defaultdict

def mutual_coverage(wordlist, threshold, concepts='concept'):
    """Compute maximal mutual coverage for all language in a wordlist.
    
    Note
    ----
    Returns all languages in a sample for which coverage is minimally as
    defined in the threshold. Coverage means the number of concepts for which
    there is a translation in both language A and language B.
    """
    def mutual(taxA, taxB):
        return len(set(
            [w for w in wordlist.get_list(
                col=taxA, 
                flat=True, 
                entry=concepts
                ) if w in wordlist.get_list(
                    col=taxB, 
                    flat=True,
                    entry=concepts)]))

    G = nx.Graph()
    for tax in wordlist.cols:
        G.add_node(tax)
    for taxA, taxB in combinations(wordlist.cols, r=2):
        coverage = mutual(taxA, taxB)
        if coverage >= threshold:
            G.add_edge(taxA, taxB, coverage=coverage)
    
    best_cliques = defaultdict(list)
    for clique in find_cliques(G):
        sums = []
        for taxA, taxB in combinations(clique, r=2):
            sums += [G.edge[taxA][taxB]['coverage']]
        if sums:
            val = int(sum(sums) / len(sums) + 0.5)
            best_cliques[val] += [sorted(clique)]
    return best_cliques
