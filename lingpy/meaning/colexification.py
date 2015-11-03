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


def _get_colexifications(wordlist, entry='ipa', concept='concept',
        family='family'):
    """
    Helper function computes colexifications for a given set of languages in a
    wordlist.
    """
    if family not in wordlist.header:
        family = 'doculect'

    taxa = wordlist.cols
    colexifications = []
    for taxon in taxa:
        
        print('Analyzing taxon {0}...'.format(taxon)) # XXX replace by log info
        
        tmp_idxs = wordlist.get_list(taxon=taxon, flat=True)
        tmp_family = wordlist[tmp_idxs[0], family]
        tmp_concepts = wordlist.get_list(taxon=taxon, flat=True,
                entry=concept)
        tmp_entries = wordlist.get_list(taxon=taxon, flat=True,
                entry=entry)
        
        # iterate over all concepts and add them to the graph
        for i,c1 in enumerate(tmp_concepts):
            for j,c2 in enumerate(tmp_concepts):
                if i < j:
                    if tmp_entries[i] == tmp_entries[j] and c1 != c2: 
                        colexifications += [(c1,c2,taxon,tmp_family,tmp_entries[i])]

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

    taxa = wordlist.cols
    statistics = {}
    for k in wordlist:
        tmp_concept = wordlist[k,concept]
        tmp_entry = wordlist[k,entry]
        tmp_family = wordlist[k,family]
        tmp_taxon = wordlist[k,'taxon']

        try:
            statistics[tmp_concept]['families'] += [tmp_family]
            statistics[tmp_concept]['doculects'] += [tmp_taxon]
            statistics[tmp_concept]['words'] += [tmp_entry]
        except KeyError:
            statistics[tmp_concept] = dict(
                    families=[tmp_family],
                    doculects=[tmp_taxon],
                    words=[tmp_entry])
    for k in statistics:
        statistics[k]['familyOcc'] = len(set(statistics[k]['families']))
        statistics[k]['doculectOcc'] = len(set(statistics[k]['doculects']))
        statistics[k]['wordOcc'] = len(statistics[k]['words'])
        statistics[k]['families'] = sorted(set(statistics[k]['families']))
        statistics[k]['doculects'] = sorted(set(statistics[k]['doculects']))

    return statistics

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

def _make_graph(colexifications, bipartite=False):
    """
    Return a graph-object from colexification data.
    """
    G = nx.Graph()

    if not bipartite:
        for c1,c2,t,f,entry in colexifications:
            try:
                G.edge[c1][c2]['families'] += [f]
                G.edge[c1][c2]['doculects'] += [t]
                G.edge[c1][c2]['words'] += [entry]
            except:
                G.add_node(c1, ntype='concept')
                G.add_node(c2, ntype='concept')
                G.add_edge(c1,c2, families=[f], doculects=[t], words=[entry])
        for a,b,d in G.edges(data=True):
            d['familyWeight'] = len(set(d['families']))
            d['wordWeight'] = len(d['words'])
            d['doculectWeight'] = len(set(d['doculects']))
            d['family'] = sorted(set(d['families']))
            d['doculects'] = sorted(set(d['doculects']))
            
    elif bipartite:
        idx = 1
        for c1,c2,t,f,entry in colexifications:
            try:
                G.edge[t+'.'+str(idx)][c1]['weight'] += 1
                G.edge[t+'.'+str(idx)][c2]['weight'] += 1
                idx += 1
            except:
                G.add_node(t+'.'+str(idx), ntype='word', entry=entry,
                        doculect=t, family=f)
                G.add_node(c1, ntype='concept')
                G.add_node(c2, ntype='concept')
                G.add_edge(t+'.'+str(idx),c1,weight=1)
                G.add_edge(t+'.'+str(idx),c2,weight=1)
                idx += 1
    
    return G

def colexification_network(wordlist, entry='ipa', concept='concept', 
        output='', filename='network', bipartite=False, **keywords):
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

    # get a list of the taxa
    taxa = wordlist.cols # XXX same as above, but we don't care in this stage

    # now, iterate over all concepts for each taxon and add the connections to
    # our network, which we now simply store as networkx graph for conveniency
    G = nx.Graph()

    colexifications = _get_colexifications(wordlist, entry, concept)
    stats = _get_statistics(wordlist, entry, concept)

    G = _make_graph(colexifications, bipartite=bipartite)

    # we should also add meta-data to the nodes in the graph
    for node,data in G.nodes(data=True):
        if data['ntype'] == 'concept':
            data.update(stats[node])

    if not output:
        return G
    
    if output == 'gml':

        for node,data in G.nodes(data=True):
            for k in data:
                if isinstance(data[k], list):
                        data[k] = '//'.join([str(x) for x in data[k]])
        for nA,nB,data in G.edges(data=True):
            for k in data:
                if isinstance(data[k], list):
                    data[k] = '//'.join([str(x) for x in data[k]])
        nx.write_gml(G, filename+'.gml')
        print("Data has been written to file {0}.".format(filename+'.gml')) # XXX

def compare_colexifications(wordlist, entry='ipa', concept='concept'):
    """
    Compare colexification patterns for a given wordlist.
    """

    taxa = wordlist.cols
    colex = {}
    
    colexifications = _get_colexifications(wordlist, entry, concept)
    colex = _get_colexifications_by_taxa(colexifications)

    return _make_matrix(taxa, colex)

def partition_colexifications(G, weight='weight'):
    """
    Use a simple partition algorithm to cluster the nodes.
    """

    # check whether keyword "weight" is actually in the graph.
    if weight != 'weight':
        for a,b,c in G.edges(data=True):
            if 'weight' not in c:
                c['weight'] = c[weight]
                
    partition = community.best_partition(G)
    converted = {}
    for node in partition:
        try:
            converted[partition[node]] += [node]
        except KeyError:
            converted[partition[node]] = [node]

    return converted, partition

def evaluate_colexifications(G, weight='word_weight', outfile=None):
    """
    Function calculates most frequent colexifications in a wordlist.

    """

    # get edges first, sorted by degree
    edges = [(a,b,c[weight]) for a,b,c in sorted(
        G.edges(data=True), key=lambda x: x[2][weight])]

    # now get the degree for the nodes
    nodes = sorted(G.degree(weight=weight).items(), key=lambda x: x[1])

    if not outfile:
        return nodes,edges

