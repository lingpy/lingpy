# *-* coding: utf-8 *-*
"""
Module provides a class for partial cognate detection, expanding the LexStat class.

"""
from __future__ import print_function, division, unicode_literals
from collections import defaultdict
from itertools import combinations_with_replacement, product, combinations
from collections import OrderedDict

import numpy as np
import networkx as nx

import lingpy
from lingpy.algorithm import clustering, extra
from lingpy.compare.lexstat import LexStat
from lingpy.util import combinations2
try:
    from lingpy.algorithm.cython import calign
except ImportError:
    from lingpy.algorithm.cython import _calign as calign

# taking functions from lexstat source code here
def _charstring(id_, char='X', cls='-'):
    return '{0}.{1}.{2}'.format(id_, char, cls)

def _get_slices(tokens, **keywords):
    """
    Function returns a slice tuple that indicates the boundaries for words with\
            morpheme annotations.
    """
    kw = dict(
            separators='_#+'+lingpy.settings.rcParams['morpheme_separator'],
            tones=lingpy.settings.rcParams['tones']
            )
    kw.update(keywords)
    out = []
    start = 0
    for i,token in enumerate(tokens):
        if token in kw['separators']:
            if tokens[i-1][0] in kw['tones']:
                start = i+1
            else:
                out += [(start, i)]
                start = i+1
                
        if token[0] in lingpy.settings.rcParams['tones']:
            out += [(start, i+1)]
            start = i+1
    if start < len(tokens):
        out += [(start, len(tokens))]
    return out

class Partial(LexStat):
    """
    Extended class for automatic detection of partial cognates.

    Parameters
    ----------
    filename : str
        The name of the file that shall be loaded.
    model : :py:class:`~lingpy.data.model.Model`
        The sound-class model that shall be used for the analysis. Defaults to
        the SCA sound-class model.
    merge_vowels : bool (default=True)
        Indicate whether consecutive vowels should be merged into single tokens or kept
        apart as separate tokens.
    transform : dict
        A dictionary that indicates how prosodic strings should be simplified
        (or generally transformed), using a simple key-value structure with the
        key referring to the original prosodic context and the value to the new
        value.  Currently, prosodic strings (see
        :py:meth:`~lingpy.sequence.sound_classes.prosodic_string`) offer 11
        different prosodic contexts. Since not all these are helpful in
        preliminary analyses for cognate detection, it is useful to merge
        some of these contexts into one. The default settings distinguish only
        5 instead of 11 available contexts, namely:

        * ``C`` for all consonants in prosodically ascending position,
        * ``c`` for all consonants in prosodically descending position,
        * ``V`` for all vowels,
        * ``T`` for all tones, and
        * ``_`` for word-breaks.

        Make sure to check also the "vowel" keyword when initialising a LexStat
        object, since the symbols you use for vowels and tones should be
        identical with the ones you define in your transform dictionary.
    vowels : str (default="VT\_")
        For scoring function creation using the
        :py:class:`~lingpy.compare.lexstat.LexStat.get_scorer` function, you
        have the possibility to use reduced scores for the matching of tones
        and vowels by modifying the "vscale" parameter, which is set to 0.5 as
        a default.  In order to make sure that vowels and tones are properly
        detected, make sure your prosodic string representation of vowels
        matches the one in this keyword. Thus, if you change the prosodic
        strings using the "transform" keyword, you also need to change the
        vowel string, to make sure that "vscale" works as wanted in the
        :py:class:`~lingpy.compare.lexstat.LexStat.get_scorer` function.
    check : bool (default=False)
        If set to **True**, the input file will first be checked for errors
        before the calculation is carried out. Errors will be written to the
        file **errors**, defaulting to ``errors.log``. See also ``apply_checks``
    apply_checks : bool (default=False)
        If set to **True**, any errors identified by `check` will be handled
        silently.
    no_bscorer: bool (default=False)
        If set to **True**, this will suppress the creation of a
        language-specific scoring function (which may become quite large and is
        additional ballast if the method "lexstat" is not used after all. If
        you use the "lexstat" method, however, this needs to be set to
        **False**.
    errors : str
        The name of the error log.

    Attributes
    ----------
    pairs : dict
        A dictionary with tuples of language names as key and indices as value, \
        pointing to unique combinations of words with the same meaning in all \
        language pairs.
    model : :py:class:`~lingpy.data.model.Model`
        The sound class model instance which serves to convert the phonetic
        data into sound classes.
    chars : list
        A list of all unique language-specific character types in the
        instantiated LexStat object. The characters in this list consist of

        * the language identifier (numeric, referenced as "langid" as a
          default, but customizable via the keyword "langid")
        * the sound class symbol for the respective IPA transcription value
        * the prosodic class value

        All values are represented in the above order as one string, separated
        by a dot. Gaps are also included in this collection. They are
        traditionally represented as "X" for the sound class and "-" for the
        prosodic string.
    rchars : list
        A list containing all unique character types across languages. In
        contrast to the chars-attribute, the "rchars" (raw chars) do not
        contain the language identifier, thus they only consist of two values,
        separated by a dot, namely, the sound class symbol, and the prosodic
        class value.
    scorer : dict
        A collection of :py:class:`~lingpy.algorithm.cython.misc.ScoreDict`
        objects, which are used to score the strings. LexStat distinguishes two
        different scoring functions:

        * rscorer: A "raw" scorer that is not language-specific and consists
          only of sound class values and prosodic string values. This scorer is
          traditionally used to carry out the first alignment in order to
          calculate the language-specific scorer. It is directly accessible as an
          attribute of the LexStat class
          (:py:class:`~lingpy.compare.lexstat.lexstat.rscorer`). The characters
          which constitute the values in this scorer are accessible via the
          "rchars" attribue of each lexstat class.
        * bscorer: The language-specific scorer. This scorer is made of unique
          language-specific characters. These are accessible via the "chars"
          attribute of each LexStat class. As the "rscorer", the "bscorer" can
          also be accessed directly as an attribute of the LexStat class 
          (:py:class:`~lingpy.compare.lexstat.lexstat.bscorer`).

    Notes
    -----
    This method automatically infers partial cognate sets from data which was
    previously morphologically segmented. 

    """

    def __init__(self, infile, **keywords):

        kw = {
                "morphemes" : "morphemes",
                "partial_cognates" : "partial_cognate_sets"
                }
        kw.update(keywords)
        lingpy.compare.lexstat.LexStat.__init__(self, infile, **kw)
        self._morphemes = kw['morphemes']
        self._partials = kw['partial_cognates']

    def _get_partial_matrices(
            self,
            concept=False,
            method='sca',
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            mode='global',
            gop=-2,
            restriction='',
            **keywords
            ):
        """
        Function creates matrices for the purpose of partial cognate detection.
        """

        # set the defaults
        kw = dict(
            defaults=False,
            external_scorer=False,  # external scoring function
            imap_mode= False,
            separators=lingpy.settings.rcParams['morpheme_separator']+'+_#'
        )
        kw.update(keywords)
        
        def function(idxA, idxB, sA, sB, **keywords):
            if method == 'lexstat':
                args = [
                        self[idxA, self._numbers][sA[0]:sA[1]],
                        self[idxB, self._numbers][sB[0]:sB[1]],
                        [self.cscorer[_charstring(
                            self[idxB, self._langid]
                            ), n]
                            for n in self[idxA, self._numbers][sA[0]:sA[1]]],
                        [self.cscorer[_charstring(
                            self[idxA, self._langid]), n]
                            for n in self[idxB, self._numbers][sB[0]:sB[1]]],
                        self[idxA, self._prostrings][sA[0]:sA[1]],
                        self[idxB, self._prostrings][sB[0]:sB[1]],
                        1,
                        scale,
                        factor,
                        self.cscorer,
                        mode,
                        restricted_chars,
                        1]
            elif method == 'sca':
                args = [
                        [n.split('.', 1)[1] for n in self[idxA,
                            self._numbers][sA[0]:sA[1]]],
                        [n.split('.', 1)[1] for n in self[idxB,
                            self._numbers][sB[0]:sB[1]]],
                        self[idxA, self._weights][sA[0]:sA[1]],
                        self[idxB, self._weights][sB[0]:sB[1]],
                        self[idxA, self._prostrings][sA[0]:sA[1]],
                        self[idxB, self._prostrings][sB[0]:sB[1]],
                        gop,
                        scale,
                        factor,
                        self.rscorer,
                        mode,
                        restricted_chars,
                        1]
            return calign.align_pair(*args)[2]
        
        concepts = [concept] if concept else sorted(self.rows)
        
        # we have two basic constraints in the algorithm:
        # a) set cognacy between morphemes in the same word to zero
        # b) set cognacy for those parts to zero which are superceded by
        # another part in all comparisons of two words
        # essentially, setting things to zero, means setting them to 1, since
        # we are dealing with distances here
        for c in concepts:
            
            indices = self.get_list(row=c, flat=True)
            matrix = []
            tracer = []
            
            # first assemble all partial parts
            trace = defaultdict(list) # stores where the stuff is in the matrix
            count = 0
            for idx in indices:
                
                # we need the slices for both words, so let's just take the
                # tokens for this time
                tokens = self[idx, self._segments]
                
                # now get the slices with the function
                slices = _get_slices(tokens, **kw)

                for i,slc in enumerate(slices):
                    tracer += [(idx, i, slc)]
                    trace[idx] += [(i, slc, count)]
                    count += 1
            
            if kw['imap_mode']:
                # now, iterate for each string pair, asses the scores, and make
                # sure, we only assign the best of those to the matrix

                matrix = [[0 for i in tracer] for j in tracer]
                # reset the self-constraints (we missed it before)


                for idxA, idxB in combinations(indices, r=2):
                    # iterate over all parts
                    scores = []
                    idxs = []
                    for i,sliceA,posA in trace[idxA]:
                        for j,sliceB,posB in trace[idxB]:
                            d = function(idxA, idxB, sliceA, sliceB)
                            scores += [d]
                            idxs += [(posA,posB)]
                    
                    visited_seqs = set([])
                    while scores:
                        min_score_index = scores.index(min(scores))
                        min_score = scores.pop(min_score_index)
                        posA, posB = idxs.pop(min_score_index)
                        if posA in visited_seqs or posB in visited_seqs:
                            matrix[posA][posB] = 1
                            matrix[posB][posA] = 1
                        else:
                            matrix[posA][posB] = min_score
                            matrix[posB][posA] = min_score
                            visited_seqs.add(posA)
                            visited_seqs.add(posB)
                for idx in indices:
                    for i,(_,sliceA,posA) in enumerate(trace[idx]):
                        for j,(_,sliceB,posB) in enumerate(trace[idx]):

                            if i < j:
                                matrix[posA][posB] = 1
                                matrix[posB][posA] = 1
            else:
                matrix = []
                for (idxA, posA, sliceA), (idxB, posB, sliceB) in combinations(tracer, r=2):
                    
                    if idxA == idxB:
                        d = 1
                    else:
                        try:
                            d = function(idxA, idxB, sliceA, sliceB)
                        except ZeroDivisionError:
                            lingpy.log.warn(
                                "Encountered Zero-Division for the comparison of "
                                "{0} and {1}".format(
                                    ''.join(self[idxA, self._tokens]),
                                    ''.join(self[idxB, self._tokens])))
                            d = 100
                    matrix += [d]
                matrix = lingpy.algorithm.misc.squareform(matrix)
            if not concept:
                yield c, tracer, matrix
            else:
                yield matrix

    def partial_cluster(
            self,
            method='sca',
            threshold=0.45,
            scale=0.5,
            factor=0.3,
            restricted_chars='_T',
            mode='overlap',
            cluster_method='infomap',
            gop=-1,
            restriction='',
            ref='',
            external_function=None,
            **keywords):
        """
        Cluster the words into partial cognate sets.

        Function for flat clustering of words into cognate sets.

        Parameters
        ----------
        method : {'sca','lexstat','edit-dist','turchin'} (default='sca')
            Select the method that shall be used for the calculation.
        cluster_method : {'upgma','single','complete', 'mcl'} (default='upgma')
            Select the cluster method. 'upgma' (:evobib:`Sokal1958`) refers to
            average linkage clustering, 'mcl' refers to the "Markov Clustering
            Algorithm" (:evobib:`Dongen2000`).
        threshold : float (default=0.3)
            Select the threshold for the cluster approach. If set to c{False},
            an automatic threshold will be calculated by calculating the
            average distance of unrelated sequences (use with care).
        scale : float (default=0.5)
            Select the scale for the gap extension penalty.
        factor : float (default=0.3)
            Select the factor for extra scores for identical prosodic segments.
        restricted_chars : str (default="T_")
            Select the restricted chars (boundary markers) in the prosodic
            strings in order to enable secondary alignment.
        mode : {'global','local','overlap','dialign'} (default='overlap')
            Select the mode for the alignment analysis.
        verbose : bool (default=False)
            Define whether verbose output should be used or not.
        gop : int (default=-2)
            If 'sca' is selected as a method, define the gap opening penalty.
        restriction : {'cv'} (default="")
            Specify the restriction for calculations using the edit-distance.
            Currently, only "cv" is supported. If *edit-dist* is selected as
            *method* and *restriction* is set to *cv*, consonant-vowel matches
            will be prohibited in the calculations and the edit distance will
            be normalized by the length of the alignment rather than the length
            of the longest sequence, as described in :evobib:`Heeringa2006`.
        inflation : {int, float} (default=2)
            Specify the inflation parameter for the use of the MCL algorithm.
        expansion : int (default=2)
            Specify the expansion parameter for the use of the MCL algorithm.
        
        """
        kw = dict(
                imap_mode = True,
                post_processing = False,
                inflation=2,
                expansion=2,
                max_steps=1000,
                add_self_loops=True,
                mcl_logs=lambda x: -np.log2((1 - x) ** 2)
                )
        kw.update(keywords)        

        # check for parameters and add clustering, in order to make sure that
        # analyses are not repeated
        if not hasattr(self, 'params'):
            self.params = {}
        self.params['partial_cluster'] = "{0}_{1}_{2:.2f}".format(
            method, cluster_method, threshold)
        self._stamp += '# Partial Cluster: ' + self.params['partial_cluster']

        matrices = self._get_partial_matrices(method=method, scale=scale,
                factor=factor, restricted_chars=restricted_chars, mode=mode,
                gop=gop, imap_mode=kw['imap_mode'])
        k = 0
        C = defaultdict(list) # stores the pcogids
        G = {} # stores the graphs
        with lingpy.util.ProgressBar('PARTIAL SEQUENCE CLUSTERING',
                len(self.rows)) as progress:
            for concept, trace, matrix in matrices:
                progress.update()
                lingpy.log.info('Analyzing concept {0}...'.format(concept))
                if external_function:
                    c = external_function(threshold, matrix,
                            taxa=list(range(len(matrix))), revert=True)
                elif cluster_method == 'infomap':
                    c = extra.infomap_clustering(threshold,
                            matrix, taxa=list(range(len(matrix))), 
                            revert=True)
                elif cluster_method == 'mcl':
                    c = clustering.mcl(threshold, matrix, 
                            taxa = list(range(len(matrix))),
                            max_steps=kw['max_steps'],
                            inflation=kw['inflation'],
                            expansion=kw['expansion'],
                            add_self_loops=kw['add_self_loops'],
                            logs=kw['mcl_logs'],
                            revert=True)
                elif cluster_method in ['upgma', 'single', 'complete', 'ward']:
                    c = clustering.flat_cluster(cluster_method,
                            threshold, matrix,
                            revert=True)
                else:
                    raise ValueError("No suitable cluster method specified.")
                
                for i,(idx,pos,slc) in enumerate(trace):
                    C[idx] += [c[i] + k]
                if kw['post_processing']:
                    _g = nx.Graph()
                    for i,(idx,pos,slc) in enumerate(trace):
                        _g.add_node((i,idx,pos))
                    remove_edges = []
                    for (i, n1), (j, n2) in combinations2(enumerate(_g.nodes())):
                        if C[n1[1]][n1[2]] == C[n2[1]][n2[2]]:
                            _g.add_edge(n1, n2)
                            if n1[1] == n2[1]:
                                # get scores for n1 and n2 with all the rest in
                                # the matrix to decide for one
                                sn1, sn2 = 0, 0
                                for i,row in enumerate(matrix):
                                    sn1 += matrix[i][n1[0]]
                                    sn2 += matrix[i][n2[0]]
                                sn1 = sn1 / len(matrix)
                                sn2 = sn2 / len(matrix)
                                if sn1 <= sn2:
                                    remove_edges += [n2]
                                else:
                                            remove_edges += [n1]
                    for node in remove_edges:
                        for edge in sorted(_g.edge[node]):
                            _g.remove_edge(node, edge)

                    for i,coms in enumerate(nx.connected_components(_g)):
                        cogid = i + 1 + k
                        for j,idx,pos in coms:
                            C[idx][pos] = cogid
                    
                    G[concept] = _g

                k += max(c.values())
        self.add_entries(ref or self._partials, C, lambda x: x)
        self.graphs = G

    def add_cognate_ids(self, source, target, idtype='strict', override=False):
        """
        Compute normal cognate identifiers from partial cognate sets.

        Parameters
        ----------
        source: str
            Name of the source column in your wordlist file.
        target : str
            Name of the target column in your wordlist file.
        idtype : str (default="strict")
            Select between "strict" and "loose".
        override: bool (default=False)
            Specify whether you want to override existing columns.
        
        Notes
        -----
        While the computation of strict cognate IDs from partial cognate IDs is
        straightforward and just judges those words as cognate which are
        identical in all their parts, the computation of loose cognate IDs
        constructs a network between all words, draws lines between all words
        that share a common morpheme, and judges all connected components in this
        network as cognate.
        """
        if idtype == 'strict':
            
            tmp = defaultdict(list)
            for k in self._data:
                tmp[tuple(self[k, source])] += [k]
            idx = 1
            D = {}
            for vals in tmp.values():
                for k in vals:
                    D[k] = idx
                idx += 1
            self.add_entries(target, D, lambda x: x, override=override)
        elif idtype == 'loose':

            D = {}
            idx = 1
            for c in self.rows:
                idxs = self.get_list(row=c, flat=True)
                srcs = [self[k, source] for k in idxs]

                # get connected components
                g = nx.Graph()
                g.add_nodes_from(idxs)
                for (i, cogsA), (j, cogsB) in combinations2(zip(idxs, srcs)):
                     if [x for x in cogsA if x in cogsB]:
                         g.add_edge(i, j)
                for i,comps in enumerate(nx.connected_components(g)):
                    for comp in comps:
                        D[comp] = idx + i
                idx += (i+1)
            self.add_entries(target, D, lambda x: x, override=override)
        else:
            raise ValueError("The value you selected is not available.")
