"""
Module offers methods to handle colexification patterns in wordlist objects.
"""
try:
    import networkx as nx
except ImportError:
    print ("Warning") # add log-specification later XXX

try:
    import community
except ImportError:
    print("Warning") # add log-specification later XXX


def _get_colexifications(wordlist, entry='ipa', concept='concept'):
    """
    Helper function computes colexifications for a given set of languages in a
    wordlist.
    """

    taxa = wordlist.cols
    colexifications = []
    for taxon in taxa:
        
        print('Analyzing taxon {0}...'.format(taxon)) # XXX replace by log info

        tmp_concepts = wordlist.get_list(taxon=taxon, flat=True,
                entry=concept)
        tmp_entries = wordlist.get_list(taxon=taxon, flat=True,
                entry=entry)
        
        # iterate over all concepts and add them to the graph
        for i,c1 in enumerate(tmp_concepts):
            for j,c2 in enumerate(tmp_concepts):
                if i < j:
                    if tmp_entries[i] == tmp_entries[j] and c1 != c2: 
                        colexifications += [(c1,c2,taxon,tmp_entries[i])]

    return colexifications

def _get_colexifications_by_taxa(colexifications):
    
    colex = {}
    for c1,c2,t,entry in colexifications:
        try:
            colex[t] += [(c1,c2),(c2,c1)]
        except KeyError:
            colex[t] = [(c1,c2),(c2,c1)]
    for t in colex:
        colex[t] = set(colex[t])
    return colex

def _make_matrix(taxa, colex):
    """
    Take colexification data and use it to create a distance matrix.

    Note
    ----
    "colex" is a dictionary with taxon names as keys and colexification data in
    form of tuples of concepts, not necessarily ordered, in both directions, as
    values.
    """
    
    # calculate the matrix
    matrix = [[0 for i in range(len(colex))] for j in range(len(colex))]
    for i,t1 in enumerate(taxa):
        for j,t2 in enumerate(taxa):
            if i < j:
                intersection = colex[t1].intersection(colex[t2])
                union = colex[t1].union(colex[t2])

                dist = 1 - len(intersection) / len(union)
                
                matrix[i][j] = dist
                matrix[j][i] = dist
    return matrix

def _make_graph(colexifications):
    """
    Return a graph-object from colexification data.
    """
    G = nx.Graph()
    for c1,c2,t,entry in colexifications:
        try:
            G.edge[c1][c2]['weight'] += 1
            G.edge[c1][c2]['occs'] += [t]
            G.edge[c1][c2]['entries'] += [entry]
        except:
            G.add_edge(c1,c2, weight=1, occs=[t], entries=[entry])
    
    return G

def colexification_network(wordlist, entry='ipa', concept='concept', 
        output=None, **keywords):
    """
    Calculate a colexification networkw from a given wordlist object.

    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        The wordlist object containing the data.

    entry : str (default="ipa")
        The reference point for the language entry. We use "ipa" as a default.
    concept : str (default="concept")
        The reference point for the name of the row containing the concepts. We
        use "concept" as a default.
    """

    # get a list of the taxa
    taxa = wordlist.cols # XXX same as above, but we don't care in this stage

    # now, iterate over all concepts for each taxon and add the connections to
    # our network, which we now simply store as networkx graph for conveniency
    G = nx.Graph()

    colexifications = _get_colexifications(wordlist, entry, concept)
    G = _make_graph(colexifications)

    if not output:
        return G

def compare_colexifications(wordlist, entry='ipa', concept='concept'):
    """
    Compare colexification patterns for a given wordlist.
    """

    taxa = wordlist.cols
    colex = {}
    
    colexifications = _get_colexifications(wordlist, entry, concept)
    colex = _get_colexifications_by_taxa(colexifications)

    return _make_matrix(taxa, colex)

def partition_colexifications(G):

    partition = community.best_partition(G)
    converted = {}
    for node in partition:
        try:
            converted[partition[node]] += [node]
        except KeyError:
            converted[partition[node]] = [node]

    return converted, partition


