# *-* coding: utf-8 *-*
"""
Module for the derivation of sound class models.

The module provides functions for the customized compilation of sound-class models.
All models are defined in simple text files. In order to guarantee their quick
access when loading the library, the models are compiled and stored in binary
files.
"""
from __future__ import division, print_function, unicode_literals
import os

from lingpy.algorithm import misc
from lingpy.convert.strings import scorer2str
from lingpy.read.phylip import read_scorer
from lingpy import cache
from lingpy import util
from lingpy import log

import networkx as nx

def _read(filename, normalize=None):
    res = {}
    for line in util.read_text_file(filename, normalize=normalize, lines=True):
        k, v = line.split(' : ')
        res[k] = v.split(', ')
    return res


def _import_sound_classes(filename):
    """
    Function imports individually defined sound classes from a text file and
    creates a replacement dictionary from these sound classes.
    """
    sc_repl_dict = {}
    errors = []
    for key, values in _read(filename, normalize='NFC').items():
        for value in values:
            log.info('%s' % ((value, key),))
            if value in sc_repl_dict and sc_repl_dict[value] != key:
                errors += [value]
            sc_repl_dict[value] = key
    if errors:
        raise ValueError("Values {0} in file {1} are multiply defined!".format(
            ' // '.join(sorted(set(errors))), filename))
    return sc_repl_dict


def _import_score_tree(filename):
    """
    Function imports score trees for a given range of sound classes and
    converts them into a graph.
    """
    graph = nx.DiGraph()
    for key, values in _read(filename).items():
        graph.add_node(key, val=values[0])
        for value in values[1:]:
            if value != '-':
                node, weight = value.split(':')
                graph.add_edge(key, node, weight=int(weight))
    return graph


def _fop(graph, start, end, path=[]):
    """
    Function returns all paths (_fop=find_all_paths) which connect to nodes in a network.
    """
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph.node:
        return []
    paths = []
    for node in graph.edge[start].keys():
        if node not in path:
            newpaths = _fop(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


def _find_dir_path(graph, start, end):
    """
    Function finds the path connecting two nodes in a directed graph under the
    condition that the two nodes are connected either directly or by a common
    ancestor node.
    """

    # first possibility: there is a direct path between the two nodes
    # if nx.shortest_path(graph,start,end) != False:
    try:
        check = nx.shortest_path(graph, start, end)
    except:
        check = False

    if check == False:

        # return nx.shortest_path(graph,start,end)
        # else:
        # except:
        # second possibility: there is a direct path between the two nodes, but
        # it starts from the other node
        # if nx.shortest_path(graph,end,start) != False:
        try:
            check = nx.shortest_path(graph, end, start)
        except:
            check = False
            # return nx.shortest_path(graph,end,start)
        # third possibility: there is no direct path between the nodes in
        # neither direction, but there is a path in an undirected graph
        if check == False:
            if _fop(graph.to_undirected(), start, end) != []:
                # here, we simply check, whether with in all paths connecting the
                # two nodes there is a node which directly connects to both nodes
                # (i.e. which is the ancestor of both nodes). If this is the case,
                # the respective shortest path is what we are looking for.
                paths = _fop(graph.to_undirected(), start, end)
                current_path_length = max([len(path) for path in paths])
                shortest_paths = nx.shortest_path(graph)
                current_path = []
                for path in paths:
                    for node in path[1:-1]:
                        if start in shortest_paths[node].keys() \
                                and end in shortest_paths[node].keys():
                            if len(path) <= current_path_length:
                                current_path_length = len(path)
                                current_path = path
                                break
                if current_path != []:
                    return current_path
                else:
                    return False
            # fourth condition: there is no path connecting the nodes at all
            else:
                return False
        else:
            return check
    else:
        return check


def _get_path_length(graph, path):
    """
    Function returns the length of a path in a weighted graph.
    """

    if path == False:
        return False
    edges = zip(path[:-1], path[1:])
    counter = 0
    for node1, node2 in edges:
        counter += graph.to_undirected()[node1][node2]['weight']
    return counter


def _make_scoring_dictionary(graph):
    """
    Function creates a scoring dictionary for individually defined sound
    classes and individually created scoring trees by counting the path length
    connecting all nodes and assigning different start weights for vowels and
    consonants.
    """
    # the scoring dictionary which will be returned by the function
    score_dict = {}

    def get_starting_value(graph, node1, node2, max_, default):
        distance = _get_path_length(
            graph, _find_dir_path(graph, node1, node2))
        # make sure that the distance doesn't exceed the default value.
        if distance is False or distance > max_:
            return default
        return max_ + default - distance

    # iterate over all nodes in the previously created graph of sound class
    # transitions
    for node1 in graph.nodes():
        for node2 in graph.nodes():
            # check, whether the key has already been created
            try:
                score_dict[(node1, node2)]
            # if not, create the key
            except KeyError:
                # if the nodes are the same, assign them the values for
                # vowel-vowel or consonant-consonant identity
                # these values might be made changeable in later versions
                if node1 == node2:
                    # for vowels and glides, the same starting value is assumed
                    if graph.node[node1]['val'] in ['v', 'g']:
                        value = 5
                    # make sure, that tones do not score
                    elif graph.node[node1]['val'] == 't':
                        value = 2
                    else:
                        value = 10
                # if the nodes are different, see, if there is a connection
                # between them defined in the directed network
                else:
                    # treat vowel-vowel and consonant-consonant matches
                    # differently
                    if graph.node[node1]['val'] == graph.node[node2]['val']:

                        # for vowels and glides, the starting value to subtract the
                        # weighted pathlength from is the vowel-vowel-identity
                        # score
                        if graph.node[node1]['val'] in ['v', 'g']:
                            # make sure that the distance doesn't exceed the
                            # default value for vowel-vowel matches, which
                            # should be zero, if there is no connection in the
                            # path defined
                            value = get_starting_value(graph, node1, node2, 5, 0)

                        # for consonants, the starting value is the
                        # consonant-consonant score
                        elif graph.node[node1]['val'] == 'c':
                            # make sure that the minimum value of C-C-matches is zero
                            value = get_starting_value(graph, node1, node2, 10, 0)
                        else:
                            # make sure that tone-tone classes score with zero
                            value = 1
                    # for vowel-consonant, vowel-glide and glide-consonant
                    # matches, the starting value is the vowel-vowel score (may
                    # also be changed in later versions)
                    else:
                        choices = [graph.node[node1]['val'], graph.node[node2]['val']]

                        # make sure to exclude tones from all matchings in
                        # order to force the algorithm to align tones with
                        # tones or gaps and with nothing else
                        if 't' in choices:
                            value = -20
                        # matches of glides with different classes
                        elif 'g' in choices:
                            # glides and vowels or glides and consonants
                            if 'v' in choices or 'c' in choices:
                                value = get_starting_value(graph, node1, node2, 10, -5)
                            else:
                                raise ValueError
                        else:
                            value = get_starting_value(graph, node1, node2, 15, -10)

                score_dict[(node1, node2)] = value
                score_dict[(node2, node1)] = value

    # add the characters for gaps in the multiple alignment process
    # note that gaps and gaps should be scored by zero according to Feng &
    # Doolittle. so far I have scored them as -1, and scoring gaps as zero made
    # the alignments getting worse, probably because most tests have been based
    # on profiles. we probably need a very good gap score.
    for node in graph.nodes():
        # missing data
        score_dict[(node, '0')] = 0
        score_dict[('0', node)] = 0

        # swaps
        score_dict[(node, '+')] = -100
        score_dict[('+', node)] = -100

        # specific values
        if graph.node[node]['val'] == 'v':
            score_dict[(node, 'X')] = 0
            score_dict[('X', node)] = 0
        elif graph.node[node]['val'] == 'g':
            score_dict[(node, 'X')] = 0
            score_dict[('X', node)] = 0
        else:
            score_dict[(node, 'X')] = 0
            score_dict[('X', node)] = 0

    score_dict[('X', '+')] = -5
    score_dict[('+', 'X')] = -5
    score_dict[('+', '+')] = 0
    score_dict[('0', '0')] = 0
    score_dict[('0', 'X')] = 0
    score_dict[('X', '0')] = 0

    # define the gaps
    score_dict[('X', 'X')] = 0
    return score_dict


def _export_score_dict(score_dict):
    """
    Function exports a scoring dictionary to a csv-file.

    @todo: This function can be better ported to another file.
    """
    letters = list(set([key[0] for key in score_dict.keys()]))
    rows = [['+'] + letters]
    for l1 in letters:
        rows.append([l1] + [str(score_dict[(l1, l2)]) for l2 in letters])
    util.write_text_file('score_dict.csv', '\n'.join('\t'.join(row) for row in rows))


def compile_model(model, path=None):
    """
    Function compiles customized sound-class models.

    Parameters
    ----------

    model : str
        A string indicating the name of the model which shall be created.

    path : str
        A string indication the path where the model-folder is stored.

    Notes
    -----
    A model is defined by a folder placed in :file:`data/models` directory of
    the LingPy package. The name of the folder reflects the name of the model.
    It contains three files: the file :file:`converter`, the file :file:`INFO`,
    and the optional file :file:`scorer`. The format requirements for these
    files are as follows:

    :file:`INFO`
        The ``INFO``-file serves as a reference for a given sound-class model.
        It can contain arbitrary information (and also be empty). If one wants
        to define specific characteristics, like the ``source``, the
        ``compiler``, the ``date``, or a ``description`` of a given model,
        this can be done by employing a key-value structure in which the key is
        preceded by an ``@`` and followed by a colon and the value is written
        right next to the key in the same line, e.g.::
            
            @source: Dolgopolsky (1986)

        This information will then be read from the ``INFO`` file and rendered
        when printing the model to screen with help of the :py:func:`print`
        function.

    :file:`converter`
        The ``converter`` file contains all sound classes which are matched
        with their respective sound values. Each line is reserved for one
        class, precede by the key (preferably an ASCII-letter) representing the
        class::

            B : ɸ, β, f, p͡f, p͜f, ƀ
            E : ɛ, æ, ɜ, ɐ, ʌ, e, ᴇ, ə, ɘ, ɤ, è, é, ē, ě, ê, ɚ
            D : θ, ð, ŧ, þ, đ
            G : x, ɣ, χ
            ...
    
    :file:`matrix`
        A scoring matrix indicating the alignment scores of all sound-class
        characters defined by the model. The scoring is structured as a simple
        tab-delimited text file. The first cell contains the character names,
        the following cells contain the scores in redundant form (with both
        triangles being filled)::

            B  10.0 -10.0   5.0 ...
            E -10.0   5.0 -10.0 ...
            F   5.0 -10.0  10.0 ...
            ...

    :file:`scorer`
        The ``scorer`` file (which is optional) contains the graph of
        class-transitions which is used for the calculation of the scoring
        dictionary. Each class is listed in a separate line, followed by the
        symbols ``v``,``c``, or ``t`` (indicating whether the class
        represents vowels, consonants, or tones), and by the classes it is
        directly connected to. The strength of this connection is indicated by
        digits (the smaller the value, the shorter the path between the
        classes)::

            A : v, E:1, O:1
            C : c, S:2
            B : c, W:2
            E : v, A:1, I:1
            D : c, S:2
            ...
        
        The information in such a file is automatically converted into a
        scoring dictionary (see :evobib:`List2012b` for details).

    Based on the information provided by the files, a dictionary for the
    conversion of IPA-characters to sound classes and a scoring dictionary are
    created and stored as a binary.  The model can be loaded with help of the
    :py:class:`~lingpy.data.model.Model` class and used in the various classes
    and functions provided by the library.
    
    See also
    --------
    lingpy.data.model.Model
    compile_dvt

    """
    log.info("Compiling model <" + model + ">...")
    # get the path to the models
    new_path = lambda *cmps: os.path.join(path or util.data_path('models'), model, *cmps)

    log.debug("Model-Path: %s" % new_path)

    # load the sound classes
    sound_classes = _import_sound_classes(new_path('converter'))

    # dump the data
    cache.dump(sound_classes, model + '.converter')
    log.info("... successfully created the converter.")

    # try to load the scoring function or the score tree
    scorer = False

    if os.path.isfile(new_path('matrix')):
        scorer = read_scorer(new_path('matrix'))
    elif os.path.isfile(new_path('scorer')):
        score_tree = _import_score_tree(new_path('scorer'))

        # calculate the scoring dictionary
        score_dict = _make_scoring_dictionary(score_tree)

        # make score_dict a ScoreDict instance
        chars = sorted(set([s[0] for s in score_dict.keys()]))
        matrix = [[0 for i in range(len(chars))] for j in
                  range(len(chars))]
        for (i, charA), (j, charB) in util.multicombinations2(enumerate(chars)):
            if i < j:
                matrix[i][j] = score_dict.get((charA, charB), -100)
                matrix[j][i] = score_dict.get((charB, charA), -100)
            elif i == j:
                matrix[i][j] = score_dict[charA, charB]

        scorer = misc.ScoreDict(chars, matrix)
        util.write_text_file(new_path('matrix'), scorer2str(scorer))

    if scorer:
        cache.dump(scorer, model + '.scorer')
        log.info("... successfully created the scorer.")
    else:
        log.info("... no scoring dictionary defined.")

    log.info("Model <" + model + "> was compiled successfully.")


def compile_dvt(path=''):
    """
    Function compiles diacritics, vowels, and tones.

    Notes
    -----
    Diacritics, vowels, and tones are defined in the :file:`data/models/dv/` directory
    of the LingPy package and automatically loaded when loading the LingPy
    library. The values are defined as the constants
    :py:obj:`rcParams['vowels']`, :py:obj:`rcParams['diacritics']`, and
    :py:obj:`rcParams['tones']`. Their core purpose is to guide the
    tokenization of IPA strings (cf.
    :py:func:`~lingpy.sequence.sound_classes.ipa2tokens`). In order to change the
    variables, one simply has to change the text files :file:`diacritics`,
    :file:`tones`, and
    :file:`vowels` in the :file:`data/models/dv` directory. The structure of
    these files is fairly simple: Each line contains a vowel or a diacritic
    character, whereas diacritics are preceded by a dash.
    
    See also
    --------
    lingpy.data.model.Model
    lingpy.data.derive.compile_model
    """
    log.info("Compiling diacritics and vowels...")

    # get the path to the models
    if not path:
        file_path = util.data_path('models', 'dvt')
    elif path in ['evolaemp', 'el']:
        file_path = util.data_path('models', 'dvt_el')
    else:
        file_path = path

    def _read_string(name):
        # normalize stuff
        # TODO: this is potentially dangerous and it is important to decide whether
        # TODO: switching to NFD might not be a better choice
        return util.read_text_file(
            os.path.join(file_path, name), normalize='NFC').replace('\n', '')

    diacritics = _read_string('diacritics').replace('-', '')
    vowels = ''.join([v for v in _read_string('vowels') if v not in diacritics])
    tones = _read_string('tones')

    dvt = (diacritics, vowels, tones)

    if path in ['evolaemp', 'el']:
        cache.dump(dvt, 'dvt_el')
    else:
        cache.dump(dvt, 'dvt')

    log.info("Diacritics and sound classes were successfully compiled.")
