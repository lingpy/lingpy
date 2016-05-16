"""
Module offers methods to handle colexification patterns in wordlist objects.
"""
from collections import defaultdict

import networkx as nx

from lingpy import log
from lingpy.util import combinations2, dotjoin, join

def _get_colexifications(wordlist, entry='ipa', concept='concept', family='family'):
    """
    Helper function computes colexifications for a given set of languages in a
    wordlist.
    """
    if family not in wordlist.header:
        family = 'doculect'

    taxa = wordlist.cols
    colexifications = []
    for taxon in taxa:
        log.info('Analyzing taxon {0}...'.format(taxon))

        tmp_idxs = wordlist.get_list(taxon=taxon, flat=True)
        tmp_family = wordlist[tmp_idxs[0], family]
        tmp_concepts = wordlist.get_list(taxon=taxon, flat=True, entry=concept)
        tmp_entries = wordlist.get_list(taxon=taxon, flat=True, entry=entry)

        # iterate over all concepts and add them to the graph
        for (i, c1), (j, c2) in combinations2(enumerate(tmp_concepts)):
            if tmp_entries[i] == tmp_entries[j] and c1 != c2:
                colexifications += [(c1, c2, taxon, tmp_family, tmp_entries[i])]

    return colexifications


def _get_statistics(wordlist, entry='ipa', concept='concept', family='family'):
    """
    Assemble data on occurrence frequency of each concept in the data.

    Returns
    -------
    statistics : dict
        A dictionary with concepts as primary keys and dictionaries as values.
        The dictionary-values themselves yield various statistics, including
        the number of language families in which each concept occurs, the
        number of languages, and the respective words themselves.
    """
    # check whether "family" is in the wordlist.header
    if family not in wordlist.header:
        family = 'doculect'

    statistics = defaultdict(lambda: defaultdict(list))
    for k in wordlist:
        tmp_concept = wordlist[k, concept]
        statistics[tmp_concept]['families'].append(wordlist[k, family])
        statistics[tmp_concept]['doculects'].append(wordlist[k, 'taxon'])
        statistics[tmp_concept]['words'].append(wordlist[k, entry])

    for k in statistics:
        statistics[k]['familyOcc'] = len(set(statistics[k]['families']))
        statistics[k]['doculectOcc'] = len(set(statistics[k]['doculects']))
        statistics[k]['wordOcc'] = len(statistics[k]['words'])
        statistics[k]['families'] = sorted(set(statistics[k]['families']))
        statistics[k]['doculects'] = sorted(set(statistics[k]['doculects']))

    return statistics


def _get_colexifications_by_taxa(colexifications):
    colex = defaultdict(set)
    for c1, c2, t, fam, entry in colexifications:
        colex[t].add((c1, c2))
        colex[t].add((c2, c1))
    return colex


def _make_matrix(taxa, colex):
    """
    Take colexification data and use it to create a distance matrix.

    Notes
    -----
    "colex" is a dictionary with taxon names as keys and colexification data in
    form of tuples of concepts, not necessarily ordered, in both directions, as
    values.
    """
    # calculate the matrix
    matrix = [[0 for i in range(len(colex))] for j in range(len(colex))]
    for (i, t1), (j, t2) in combinations2(enumerate(taxa)):
        intersection = colex[t1].intersection(colex[t2])
        union = colex[t1].union(colex[t2])
        matrix[i][j] = matrix[j][i] = 1 - len(intersection) / len(union)
    return matrix


def _make_graph(colexifications, bipartite=False):
    """
    Return a graph-object from colexification data.
    """
    G = nx.Graph()

    if not bipartite:
        for c1, c2, t, f, entry in colexifications:
            try:
                G.edge[c1][c2]['families'] += [f]
                G.edge[c1][c2]['doculects'] += [t]
                G.edge[c1][c2]['words'] += [entry]
            except:
                G.add_node(c1, ntype='concept')
                G.add_node(c2, ntype='concept')
                G.add_edge(c1, c2, families=[f], doculects=[t], words=[entry])
        for a, b, d in G.edges(data=True):
            d['familyWeight'] = len(set(d['families']))
            d['wordWeight'] = len(d['words'])
            d['doculectWeight'] = len(set(d['doculects']))
            d['family'] = sorted(set(d['families']))
            d['doculects'] = sorted(set(d['doculects']))
    elif bipartite:
        for idx, (c1, c2, t, f, entry) in enumerate(colexifications):
            nindex = dotjoin(t, idx + 1)
            try:
                G.edge[nindex][c1]['weight'] += 1
                G.edge[nindex][c2]['weight'] += 1
            except KeyError:
                G.add_node(nindex, ntype='word', entry=entry, doculect=t, family=f)
                G.add_node(c1, ntype='concept')
                G.add_node(c2, ntype='concept')
                G.add_edge(nindex, c1, weight=1)
                G.add_edge(nindex, c2, weight=1)
    return G


def colexification_network(
        wordlist,
        entry='ipa',
        concept='concept',
        output='',
        filename='network',
        bipartite=False,
        **keywords):
    """
    Calculate a colexification network from a given wordlist object.

    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        The wordlist object containing the data.

    entry : str (default="ipa")
        The reference point for the language entry. We use "ipa" as a default.
    concept : str (default="concept")
        The reference point for the name of the row containing the concepts. We
        use "concept" as a default.
    output: str (default='')
        If output is set to "gml", the resulting network will be written to a
        textfile in GML format.

    Returns
    -------
    G : networkx.Graph
        A networkx.Graph object.

    """
    # now, iterate over all concepts for each taxon and add the connections to
    # our network, which we now simply store as networkx graph for conveniency
    colexifications = _get_colexifications(wordlist, entry, concept)
    stats = _get_statistics(wordlist, entry, concept)

    G = _make_graph(colexifications, bipartite=bipartite)

    # we should also add meta-data to the nodes in the graph
    for node, data in G.nodes(data=True):
        if data['ntype'] == 'concept':
            data.update(stats[node])

    if not output:
        return G

    def stringify_data(data):
        for k in data:
            if isinstance(data[k], list):
                data[k] = join('//', *data[k])

    if output == 'gml':
        for node, data in G.nodes(data=True):
            stringify_data(data)
        for nA, nB, data in G.edges(data=True):
            stringify_data(data)
        nx.write_gml(G, filename + '.gml')
        log.file_written(filename + '.gml')


def compare_colexifications(wordlist, entry='ipa', concept='concept'):
    """
    Compare colexification patterns for a given wordlist.
    """
    colexifications = _get_colexifications(wordlist, entry, concept)
    return _make_matrix(wordlist.cols, _get_colexifications_by_taxa(colexifications))

def evaluate_colexifications(G, weight='wordWeight', outfile=None):
    """
    Function calculates most frequent colexifications in a wordlist.

    """

    # get edges first, sorted by degree
    edges = [(a, b, c[weight]) for a, b, c in sorted(
        G.edges(data=True), key=lambda x: x[2][weight])]

    # now get the degree for the nodes
    nodes = sorted(dict(G.degree(weight=weight)).items(), key=lambda x: x[1])

    if not outfile:
        return nodes, edges
