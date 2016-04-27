# *-* coding: utf-8 *-*
"""
Phylogeny-based detection of borrowings in lexicostatistical wordlists.
"""
from __future__ import print_function, division, unicode_literals
import os
import itertools

from six import text_type

import numpy as np

from lingpy.compare._phylogeny.utils import get_acs
from lingpy.compare._phylogeny._settings import rcParams
from lingpy.align.multiple import Multiple
from lingpy.convert.plot import plot_tree
from lingpy import compat
from lingpy import util
from lingpy import log

import networkx as nx

# mpl is only used for specific plots, we can therefor make a safe import
try:  # pragma: no cover
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except:
    mpl, plt = None, None
    log.missing_module('matplotlib')

# import the geoplot module
try:  # pragma: no cover
    import mpl_toolkits.basemap as bmp
except ImportError:
    bmp = None
    log.missing_module('basemap')

try:
    import scipy as sp
except ImportError:
    sp = None
    log.missing_module('scipy')

from ._phylogeny.polygon import getConvexHull

# lingpy imports
from ..thirdparty import cogent as cg
from ..convert.graph import gls2gml, radial_layout
from ..basic import Wordlist
from ..read.csv import csv2dict, csv2list


def get_gls(
        paps,
        taxa,
        tree,
        gpl=1,
        weights=(1, 1),
        push_gains=True,
        missing_data=0):
    """
    Calculate a gain-loss scenario.

    Parameters
    ----------
    paps : list
        A list containing the presence-absence patterns for all leaves of the
        reference tree. Presence is indicated by 1, and absence by 0. Missing
        characters are indicated by -1.
    taxa : list
        The list of taxa (leaves of the tree).
    tree : str
        A tree in Newick-format. Taxon names should (of course) be identical
        with the names in the list of taxa.
    gpl : int
        Gains per lineage. Specify the maximal amount of gains per lineage. One
        lineage is hereby defined as one path in the tree. If set to 0, only
        one gain per lineage is allowed, if set to 1, one additional gain is
        allowed, and so on. Use with care, since this will lead to larger
        computation costs (more possibilities have to be taken care of) and can
        also be quite unrealistic.
    weights : tuple (default=(1,1))
        Specify the weights for gains and losses. Setting this parameter to
        (2,1) will penalize gain events with 2 and loss events with 1.
    push_gains : bool (default=True)
        Determine whether of a set of equally parsimonious patterns those
        should be retained that show gains closer to the leaves of the tree or
        not.
    missing_data : int (default=0)
        Determine how missing data should be represented. If set to 0
        (default), missing data will be treated in the same way as absence
        character states. If you want missing data to be accounted for in the
        algorithm, set this parameter to -1.

    Notes
    -----
    This is an enhanced version of the older approach to parsimony-based
    gain-loss mapping. The algorithm is much faster than the previous one and
    also written much clearer as to the code. In most tests I ran so far, it
    also outperformed other approaches by finding more parsimonious solutions.

    """
    # this line is just to make sure we actually copy the paps and don't change
    # them unwillingly
    pap = [p for p in paps]

    for i, taxon in enumerate(taxa):
        if pap[i] == -1:
            pap[i] = missing_data

    # get dictionary for taxa with their states
    statesD = dict(zip(taxa, pap))

    states1 = [s for s in statesD if statesD[s] == 1]

    # get subtree for taxa with positive paps
    tree = tree.lowestCommonAncestor([t for t in taxa if statesD[t] == 1])

    # get the tips for the current subtree
    tips = tree.getTipNames()

    # get the root name
    root = tree.Name

    # get the nodes
    nodes = tree.getNodeNames()

    # get the distance from each taxon to the root
    distances = {}
    for node in [n for n in nodes if n != root]:

        # get distance to root
        d = len(tree.getConnectingEdges(root, node))

        # append to dictionary
        try:
            distances[d] += [node]
        except:
            distances[d] = [node]
    distances[0] = [root]

    # assign the scenarios, each scenario consists of the state of the node in
    # the tree and a dictionary with the previous events, where the node-name
    # is the key and the event (1,0,-1) is the value
    scenarios = {}
    for taxon in tips:
        scenarios[taxon] = [(statesD[taxon], {})]

    # return simple scenario if the group is single origin
    if [statesD[tip] for tip in tips].count(1) == len(tips):
        return [(root, 1)]

    # start iteration over outmost layer
    for i in sorted(distances, reverse=True):
        log.debug("Calculating Layer {0}...".format(i))

        # iterate over all nodes in the layer
        for node in distances[i]:

            # get the tree node
            tree_node = tree.getNodeMatchingName(node)

            # check whether node is tip
            if not tree_node.isTip():

                # get the nodes' children
                names = [n.Name for n in tree_node.Children]
                children = [scenarios[n] for n in names]
                log.debug("... current node {0} ({1})".format(node, names))

                # get all possible combinations
                combinations = itertools.product(*children)

                # define new nodes list (to be appended to new node
                new_nodes = []

                for combination in combinations:

                    # get stories
                    states = [node[0] for node in combination]
                    stories = [node[1] for node in combination]

                    # evaluate the states
                    s1 = states.count(1)
                    s0 = states.count(0)
                    sM = states.count(-1)

                    # get the length of the states
                    sL = len(states)

                    # get new stories
                    new_stories = {}
                    for story in stories:
                        new_stories.update(story)

                    # combine states if they evaluate to 1
                    if s1 + sM == sL:

                        # append the new combined stuff to the dictionary
                        new_nodes += [(1, new_stories)]
                        log.debug("...... 1 nodes: %s" % (new_nodes[-1],))

                    # combine states if they evaluate to 0
                    elif s0 + sM == sL:

                        # append the new combined stuff to the dictionary
                        new_nodes += [(0, new_stories)]
                        log.debug("...... 0 nodes: %s" % (new_nodes[-1],))

                    # if the both evaluate to -1, also combine them
                    elif sM == sL:
                        new_nodes += [(-1, new_stories)]

                    # append both scenarios if there's both 1 and 0
                    else:

                        # assuming origin, each node that has a 0, needs an extra
                        # origin
                        new_storiesA = new_stories.copy()
                        new_storiesB = new_stories.copy()

                        for j, state in enumerate(states):
                            if state == 1:
                                new_storiesA[names[j]] = 1
                            if state == 0:
                                new_storiesB[names[j]] = 0

                        new_nodes += [(1, new_storiesB)]
                        new_nodes += [(0, new_storiesA)]

                        log.debug("...... 01 nodes: %s" % (new_nodes[-2],))
                        log.debug("...... 01 nodes: %s" % (new_nodes[-1],))

                # evaluate the scenarios for consistency reasons,
                good_nodes = []
                minGains, minLoss = {}, {}
                for j, (state, scenario) in enumerate(new_nodes):

                    # avoid to append scenarios with more than allowed gains
                    # per lineage
                    if not (state == 1 and list(scenario.values()).count(1) > gpl):

                        # check scenarios having a loss in order to retrieve
                        # the scenario with the minimal weight, since once a
                        # loss is determined, the gains can be freely chosen
                        if state == 0:
                            gains = list(scenario.values()).count(1)
                            losses = list(scenario.values()).count(0)
                            w = gains * weights[0] + losses * weights[1]
                            try:
                                minGains[w] += [j]
                            except KeyError:
                                minGains[w] = [j]

                        # do the same for scenarios having a gain, if multiple
                        # loss-models are encountered
                        elif state == 1:
                            gains = list(scenario.values()).count(1)
                            losses = list(scenario.values()).count(0)
                            w = gains * weights[0] + losses * weights[1]
                            try:
                                minLoss[w] += [j]
                            except:
                                minLoss[w] = [j]

                # append lowest weights in gains to the list
                if minGains:
                    good_nodes += [new_nodes[idx] for idx in
                                   minGains[min(minGains)]]
                if minLoss:
                    good_nodes += [new_nodes[idx] for idx in
                                   minLoss[min(minLoss)]]

                scenarios[tree_node.Name] = good_nodes

    # select the best of all scenarios by comparing all weights
    w = len(taxa) * weights[0]

    winners = {}
    for s in scenarios[root]:
        if s[0] == 1:
            s[1][root] = s[0]

        # count the weights
        events = list(s[1].values())
        gains = events.count(1)
        losses = events.count(0)

        w = gains * weights[0] + losses * weights[1]

        try:
            winners[w] += [list(s[1].items())]
        except:
            winners[w] = [list(s[1].items())]

    # select the scenario with the hightest number of gains, if push-gains
    # option is set to true
    log.debug(
        '%s' % (
            [x for x in tree.getNodeMatchingName(root).getTipNames()
             if x not in states1],))

    return sorted(
        winners[min(winners)],
        key=lambda x: [y[1] for y in x].count(1 if push_gains else 0))[0]


class PhyBo(Wordlist):
    """
    Basic class for calculations using the TreBor method.


    Parameters
    ----------
    dataset : string
        Name of the dataset that shall be analyzed.
    tree : {None, string}
        Name of the tree file.
    paps : string (default="pap")
        Name of the column that stores the specific cognate IDs consisting
        of an arbitrary integer key and a key for the concept.
    ref : string (default="cogid")
        Name of the column that stores the general cognate ids (the
        "reference" of the analysis).
    tree_calc : {'neighbor','upgma'} (default='neighbor')
        Select the algorithm to be used for the tree calculation if no tree is
        passed with the file.
    missing : int (default=-1)
        Specify how missing data should be handled. If set to -1, missing data
        can account for both presence or absence of a cognate set in the given
        language. If set to 0, missing data is treated as absence.
    degree : int (default=100)
        The degree which is chosen for the projection of the tree layout.
    """

    # XXX generally: find a way to check whether a dataset was already loaded,
    # XXX otherwise it takes too long a time to recalculate everything

    def __init__(
        self,
        dataset,
        tree=None,
        paps='pap',
        ref='cogid',
        tree_calc='neighbor',
        output_dir=None,
        **keywords
    ):
        # TODO check for keywords, allow to load trees, etc.
        defaults = {
            'degree': 100,
            'singletons': True,
            'missing': -1,
            'change': lambda x: x ** 1.5,
            'start': 0
        }
        for k in defaults:
            keywords.setdefault(k, defaults[k])

        # check for cognates
        if 'cognates' in keywords:
            log.deprecated('cognates', 'ref')
            ref = keywords['cognates']

        # store the basename of the dataset without suffix and the identifier for paps
        self._dataset_dir = os.path.dirname(os.path.abspath(dataset))
        self._output_dir = output_dir or self._dataset_dir
        dataset_name = os.path.basename(dataset)
        self.dataset = dataset_name[:-4] if dataset_name[-4:] in ['.qlc',
                                                                  '.csv'] else dataset_name
        self._pap_string = paps

        # open csv-file of the data and store it as a word list attribute
        if os.path.isfile(dataset):
            if dataset.endswith('.csv'):
                log.deprecated('csv', 'qlc')
        else:
            raise compat.FileNotFoundError("The input file could not be found.")
        Wordlist.__init__(self, dataset, row='concept', col='doculect')

        log.info("Loaded the wordlist file.")

        # check for glossid
        if 'glid' not in self.entries:
            self._gl2id = dict(
                zip(
                    self.rows,
                    [i + 1 for i in range(len(self.rows))]
                )
            )
            self._id2gl = dict([(b, a) for a, b in self._gl2id.items()])

            f = lambda x: self._gl2id[x]

            self.add_entries(
                'glid',
                'concept',
                f
            )
        else:
            self._id2gl = dict([(int(self[k, 'glid']), self[k, 'concept']) for k in self])
            self._gl2id = dict([(self[k, 'concept'], int(self[k, 'glid'])) for k in self])

        # check for paps as attribute in the wordlist
        if paps not in self.entries:
            # define the function for conversion
            f = lambda x, y: "{0}:{1}".format(x[y[0]], x[y[1]])
            self.add_entries(
                paps,
                ref + ',glid',
                f
            )
            log.info("Created entry PAP.")

        # get the paps and the etymological dictionary
        if not hasattr(self, 'paps'):
            self.paps = self.get_paps(ref=paps, missing=keywords['missing'])
            self.etd = self.get_etymdict(ref=paps)

        log.info("Created the PAP matrix.")

        # get a list of concepts corresponding to the cogs and get the
        # singletons to be excluded from the calculation
        if not hasattr(self, 'singletons'):
            tmp = self.get_etymdict(ref=paps, entry='concept')

            # a dictionary with pap-key as key and concept as value
            self.pap2con = {}

            # list stores the singletons
            self.singletons = []

            # only calculate singletons if the option is chosen
            for key in self.paps:

                # get the names of the concepts
                concept_list = [k for k in tmp[key] if k != 0]
                concept = concept_list[0][0]
                self.pap2con[key] = concept

                # check for singletons
                if keywords['singletons']:
                    if sum([1 for p in self.paps[key] if p >= 1]) == 1:
                        self.singletons.append(key)

            # create a list of keys for faster access when iterating
            self.cogs = [k for k in self.pap2con if k not in self.singletons]
            log.info("Excluded singletons.")

        # summarize the cognate sets under their common concept

        # Load the tree, if it is not defined, assume that the treefile has the
        # same name as the dataset
        if not tree:
            if not hasattr(self, 'tree'):
                # try to load the tree first
                try:
                    self.tree = cg.LoadTree(dataset + '.tre')
                except:
                    # create it otherwise
                    self.calculate(
                        'tree',
                        ref=ref,
                        tree_calc=tree_calc,
                    )
                    log.info("Tree-file was not found, creating it now...")
                    # XXX TODO

        # if it is explicitly defined, try to load that file
        else:  # not hasattr(self,'tree'):
            self._meta['tree'] = cg.LoadTree(tree)
            log.info("Loaded the tree.")

        if isinstance(self.tree, text_type):
            self.tree = cg.LoadTree(treestring=self.tree)

        # if no good topology is given, create it automatically, using
        # the radial layout function
        gTpl = radial_layout(
            self.tree,
            filename='',
            degree=keywords['degree'],
            change=keywords['change'],
            start=keywords['start']
        )
        log.info("Calculated radial layout for the tree. ")

        self.tgraph = gTpl

        # create a couple of further attributes
        for a in ['stats', 'gls', 'dists', 'graph', 'acs']:
            if not hasattr(self, a):
                setattr(self, a, {})

    def _output_path(self, *comps, **kw):
        """A path within the output directory for the dataset.

        Note: All intermediate directories will be created unless a keyword argument
        mkdir=False is passed.

        :param comps: path components relative to the output directory.
        :return: the path.
        """
        subdir = os.path.basename(self.dataset) + '_phybo'
        path = os.path.join(self._output_dir, subdir, *comps)
        if kw.get('mkdir', True):
            dirname = os.path.dirname(path)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        return path

    def _write_file(self, name, content, log=True):
        """Write a file to the dataset-specific output directory.

        :param name: Name of the file to be written relative to the output directory.
        :param content: Content to be written.
        :param log: Flag signalling whether there should be log output or not.
        """
        util.write_text_file(self._output_path(name), content, log=log)

    def _config(self):
        if 'conf' in self._meta:
            return self._meta['conf']
        try:
            return util.jsonload(self.dataset + '.json')
        except:
            return {}

    def _get_GLS_top_down(
        self,
        pap,
        mode=1,
        missing_data=0
    ):
        """
        Infer gain-loss scenario using the method by Dagan & Martin (2007).

        """
        # check for mode
        try:
            mode = int(mode)
        except:
            raise ValueError("[i] Mode should be an integer.")

        # get the list of nodes that are not missing
        taxa, paps = [], []
        for i, taxon in enumerate(self.taxa):
            if pap[i] == -1:
                pap[i] = missing_data
            taxa += [taxon]
            paps += [pap[i]]

        # get list of taxa where pap is 1
        presents = [self.taxa[i] for i in range(len(self.taxa)) if pap[i] in (1, -1)]

        # get the subtree containing all taxa that have positive paps
        tree = self.tree.lowestCommonAncestor(
            [self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1])
        log.debug("Subtree is {0}.".format(text_type(tree)))

        # assign the basic (starting) values to the dictionary
        nodes = [t.Name for t in tree.tips()]
        log.debug("Nodes are {0}.".format(','.join(nodes)))

        if mode == 1:
            return [(tree.Name, 1)]

        # store the scenario
        scenario = []

        # make the queue
        queue = [[tree, 1]]
        while queue:

            # get tree and counter from queue
            tmp_tree, counter = queue.pop(0)

            # break if counter exceeds the mode
            if counter >= mode:
                t = tmp_tree.lowestCommonAncestor([p for p in presents if p in
                                                   tmp_tree.getTipNames()])
                if hasattr(t, 'Name'):
                    scenario += [(t.Name, 1)]
                else:
                    scenario += [(tmp_tree.Name, 0)]
            else:
                # store common names and children nodes
                commons = []

                # get tip names for checking
                tmp_names = tmp_tree.getTipNames()

                # store results for separate children
                tmp_results = []

                # now, we iterate over all children and append them to the
                # queue
                children = tmp_tree.Children
                for child in children:

                    # get lowest common ancestor
                    subtree = child.lowestCommonAncestor(
                        [p for p in presents if p in child.getTipNames()]
                    )

                    # check for tip names in subtrees
                    if hasattr(subtree, "Children"):
                        subnames = subtree.getTipNames()
                    else:
                        if child.Name in presents:
                            subnames = [child.Name]
                        else:
                            subnames = []

                    # append subnames to commons
                    commons += subnames

                    # append data to our tmp_results
                    tmp_results += [child]

                # evaluate the results
                cSet = set(commons)
                tSet = set(tmp_names)
                # check for identity and stop iteration if tips are identical
                if cSet == tSet:
                    scenario += [(tmp_tree.lowestCommonAncestor(presents).Name, 1)]
                # otherwise check for intersection and small amount of
                # differences
                # elif tmp_tree == tree.lowestCommonAncestor(commons):
                # elif cSet.issubset(tSet) and len(tSet) - len(cSet) < len(tSet) / 2:
                #    scenario += [(tmp_tree.lowestCommonAncestor(presents).Name,1)]
                # otherwise append the other results to the queue
                else:
                    add2scenario = []
                    for child in tmp_results:
                        if child.Children:
                            queue += [(child, counter + 1)]
                        else:
                            if child.Name in presents:
                                add2scenario += [child]

                    if len(add2scenario) == 2:
                        if add2scenario[0].Parent == add2scenario[1].Parent:
                            scenario += [(add2scenario[0].Parent.Name, 1)]
                        else:
                            for c in add2scenario:
                                scenario += [(c.Name, 1)]
                    else:
                        for c in add2scenario:
                            scenario += [(c.Name, 1)]

        # TODO fill the scenario with gaps
        output = []
        d = {}
        for i, taxon in enumerate(taxa):
            if paps[i] >= 1:
                d[taxon] = 1
            else:
                d[taxon] = 0

        for s in scenario:

            # append scenario to output
            output += [s]

            # get the subtree
            subtree = tree.getNodeMatchingName(s[0])

            # check whether subtree is leave
            if not subtree.Children:
                pass
            else:
                # order the internal nodes according to the number of their leaves
                ordered_nodes = sorted(
                    subtree.nontips() + [subtree], key=lambda x: len(x.tips())
                )

                # start bottom-up
                for node in ordered_nodes:

                    # get the children
                    children = node.Children

                    # store the states
                    states = []

                    # iterate over the children
                    for child in children:
                        state = d[child.Name]
                        states += [state]

                    # check for identity of states
                    if sum(states) == len(states):
                        d[node.Name] = 1
                    elif sum(states) == 0:
                        d[node.Name] = 0
                    else:
                        d[node.Name] = 1
                        for i, state in enumerate(states):
                            if state == 0:
                                output += [(children[i].Name, 0)]

        return output

    def _get_GLS(
        self,
        pap,
        mode='w',
        r=(1, 1),
        gpl=1,
        push_gains=True,
        missing_data=0
    ):
        """
        Calculate a gain-loss scenario (GLS) for a given PAP.

        Parameters
        ----------
        pap : list
            The presence/absence pattern of a given cognate-set.
        mode : str (default='w')
            The mode of the analysis, select between "w" (weights) and "r"
            (restrictio).
        r : { tuple, int } (default=(1,1))
            The weights (as binary tuple) or the restriction (an integer),
            negative restrictions mark the maximal amount of losses.
        gpl : int (default=1)
            The maximal number of gains per lineage.
        push_gains : bool (default=True)
            Indicate whether gains should be pushed to the leaves or not.
        missing_data : int (default=0)
            Indicate, how missing values should be represented in the paps. If
            set to 0, missing values will be treated as non-cognate words. If
            set to 1, missing values will be treated as potential cognates.

        """
        # make a dictionary that stores the scenario
        d = {}

        # get the list of nodes that are not missing
        taxa, paps = [], []
        for i, taxon in enumerate(self.taxa):
            if pap[i] == -1:
                pap[i] = missing_data
            taxa += [taxon]
            paps += [pap[i]]

        # get the subtree containing all taxa that have positive paps
        tree = self.tree.lowestCommonAncestor(
            [self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1])

        log.debug("Subtree is {0}.".format(text_type(tree)))

        # assign the basic (starting) values to the dictionary
        nodes = [t.Name for t in tree.tips()]
        log.debug("Nodes are {0}.".format(','.join(nodes)))

        # calculate the initial restriction value (maximal weight). This is roughly
        # spoken simply the minimal value of either all events being counted as
        # origins (case 1) or assuming origin of the character at the root and
        # counting all leaves that lost the character as single loss events (case
        # 2). In case two, the first gain of the character has to be added
        # additionally
        if mode == 'w':
            RST = min(paps.count(1) * r[0], paps.count(0) * r[1] + r[0])
        elif mode == 'r':
            RST = r

        # get maximal number of gains and losses, note that we have to include
        # missing data in a two-fold fashion here. this is probably
        # computationally not the most feasible solution. however, it is the
        # only way I can think of at the moment
        maxG = sum([1 for x in nodes if paps[taxa.index(x)] in (1, -1)])
        maxL = sum([1 for x in nodes if paps[taxa.index(x)] in (0, -1)])

        log.debug("Initial restriction threshold is {0}.".format(RST))

        # get the first state of all nodes and store the state in the
        # dictionary. note that we start from two distinct scenarios: one
        # assuming single origin at the root where all present states in the
        # leave are treated as retentions, and one assuming multiple origins,
        # where all present states in the leaves are treated as origins
        dbpaps = []
        for node in nodes:
            idx = taxa.index(node)
            if paps[idx] >= 1:
                state = 1
            else:
                state = paps[idx]
            dbpaps += [node + '/' + str(state)]

            # we append the maximally remaining possible number of gains and
            # losses to the queue dictionary and decrease it steadily once two
            # branches are merged as either loss or gain
            d[node] = [(state, [], maxG, maxL)]

        # return simple scenario, if the group is single-origin
        if sum([d[node][0][0] for node in nodes]) == len(nodes):
            return [(tree.Name, 1)]

        # order the internal nodes according to the number of their leaves
        ordered_nodes = sorted(
            tree.nontips() + [tree], key=lambda x: len(x.tips())
        )

        search_space = 0
        log.debug('The Pap to be analysed: %s' % ', '.join(dbpaps))

        # join the nodes successively
        for i, node in enumerate(ordered_nodes):
            log.debug('Node to be joined in this run: %s' % node)

            # when dealing with multifurcating trees, we have to store all
            # possible scenarios, i.e. we need to store the crossproduct of all
            # scenarios

            # get the names of the children of the nodes
            names = [x.Name for x in node.Children]

            # get the nodes with their states from the dictionary
            tmp_nodes = [d[x.Name] for x in node.Children]

            # get the cross-product of the stuff
            crossp = itertools.product(*tmp_nodes)

            newNodes = []

            # combine the histories of the items if all have the same value,
            # therefore, we first get the states in a simple list
            for cross in crossp:
                search_space += 1

                states = [x[0] for x in cross]
                stories = [x[1] for x in cross]

                # get the restriction values
                maxGains = [x[2] for x in cross]
                maxLosses = [x[3] for x in cross]

                states_1 = states.count(1)
                states_0 = states.count(0)
                states_m = states.count(-1)  # missing states

                states_len = len(states)

                # get the minimal gain and loss values
                maxGain = max(maxGains)
                maxLoss = max(maxLosses)

                # calculate the restriction value
                if mode == 'w':
                    # the following line contains some serious bug
                    # (DENKFEHLER), it works at the moment, but we should be
                    # very DAMN careful with this!
                    rst = min(maxGain * r[0] + r[1], maxLoss * r[1] + r[0])
                else:
                    rst = abs(RST)

                log.debug(
                    "... MaxG / MaxL / rst: {0} / {1} / {2}.".format(maxGain, maxLoss,
                                                                     rst))
                # log.debug("... Stories: %s" % stories)

                # combine the histories
                new_stories = []
                for x in stories:
                    new_stories += x

                # if states are identical and point to gain / presence of
                # chars, we add them directly. here we also include the number
                # of missing states: if missing states turn up, we simply treat
                # them as presence values
                if states_1 + states_m == states_len:

                    # add the histories to the queue only if their weight is
                    # less or equal to the maxWeight
                    gl = [k[1] for k in new_stories] + [1]

                    if mode == 'w':
                        weight = gl.count(1) * r[0] + gl.count(0) * r[1]
                    else:
                        if RST < 0:
                            weight = gl.count(0)
                        else:
                            weight = gl.count(1) + 1  # we need to add 1 here

                    log.debug("... state,weight: %s %s" % (gl, weight))
                    # when combining two gains, make sure that the allowed
                    # amount of gains per lineage will not be overwritten by
                    # the combination of new gains
                    gains_per_lineage = sum([1 for k in new_stories if k[1] == 1])

                    if weight <= rst and gains_per_lineage < gpl:

                        # make sure to append a smaller restriction value,
                        # since we could spare one event due to regular
                        # calculation
                        if mode == 'w':
                            newNodes.append((1, new_stories, maxGain - 1, maxLoss))
                        else:
                            # if gl.count(0) <= maxLoss:
                            newNodes.append((1, new_stories, maxGain, maxLoss))

                # if states are identical and point to absence of chars, we
                # assign them directly to the higher node. here, missing chars
                # are also included
                elif states_0 + states_m == states_len:
                    gl = [k[1] for k in new_stories]
                    log.debug("... state is 1 %s" % gl)

                    if mode == 'w':
                        weight = gl.count(1) * r[0] + gl.count(0) * r[1]
                    else:
                        if RST < 0:
                            weight = gl.count(0)
                        else:
                            weight = gl.count(1)

                    log.debug("... state,weight: %s %s" % (gl, weight))

                    if weight <= rst:
                        if mode == 'w':
                            newNodes.append((0, new_stories, maxGain, maxLoss - 1))
                        else:
                            # if gl.count(0) <= maxLoss:
                            newNodes.append((0, new_stories, maxGain, maxLoss - 1))

                # if states are both missing
                elif states_m == states_len:
                    log.debug("... all states are missing")

                    if mode == 'w':
                        newNodes.append((-1, new_stories, maxGain, maxLoss - 1))
                    else:
                        # if gl.count(0) <= maxLoss:
                        newNodes.append((-1, new_stories, maxGain, maxLoss - 1))

                # if the states are not identical, we check for both scenarios
                else:
                    log.debug("... states are different.")

                    # first scenario (tmpA) assumes origin, that is, for each node
                    # that has a 1, we add an origin to new_stories, same is
                    # for loss scenario (tmpB)
                    tmpA = [x for x in new_stories]
                    tmpB = [x for x in new_stories]
                    for c, state in enumerate(states):
                        if state == 1 or state == -1:
                            tmpA += [(names[c], 1)]
                        if state == 0 or state == -1:
                            tmpB += [(names[c], 0)]

                    # get the vectors to make it easier to retrieve the number
                    # of losses and gains
                    glA = [k[1] for k in tmpA]
                    glB = [k[1] for k in tmpB] + [1]  # don't forget adding 1 origin

                    # check the gain-loss scores
                    if mode == 'w':
                        weightA = glA.count(1) * r[0] + glA.count(0) * r[1]
                        weightB = glB.count(1) * r[0] + glB.count(0) * r[1]
                    else:
                        if RST < 0:
                            weightA = glA.count(0)
                            weightB = glB.count(0)
                        else:
                            weightA = glA.count(1)
                            weightB = glB.count(1)

                    # create the new nodes. Note that we can only reduce the
                    # number of losses here by one, but not the possible number
                    # of gains, since the last gain will also score in our
                    # calculation, but we cannot predict, whether a given gain
                    # is indeed the last one (or can we?)
                    newNodeA = (0, tmpA, maxGain, maxLoss)
                    newNodeB = (1, tmpB, maxGain, maxLoss)

                    # check for additional gains in the gain-scenario,
                    # according to the current model, we don't allow for one
                    # character to be gained twice along a branch, i.e. by an
                    # ancestor, then get lost, and than be gained anew
                    gains_per_lineage = sum([1 for k in tmpB if k[1] == 1])
                    if gains_per_lineage >= gpl:  # in [k[1] for k in tmpB]:
                        noB = True
                    else:
                        noB = False

                    if weightA <= rst:
                        newNodes += [newNodeA]

                    if weightB <= rst and not noB:
                        newNodes += [newNodeB]

                d[node.Name] = newNodes
                log.debug("... Possible scenarios for '{0}': {1}".format(node.Name, len(
                    d[node.Name])))

        # try to find the best scenario by counting the ratio of gains and losses.
        # the key idea here is to reduce the number of possible scenarios according
        # to a given criterion. We choose the criterion of minimal changes as a
        # first criterion to reduce the possibilities, i.e. we weight both gains
        # and losses by 1 and select only those scenarios where gains and losses
        # sum up to a minimal number of gains and losses. This pre-selection of
        # scenarios can be further reduced by weighting gains and losses
        # differently. So in a second stage we choose only those scenarios where
        # there is a minimal amount of gains.

        # convert the specific format of the d[tree.Name] to simple format
        gls_list = []
        for first, last, mg, ml in d[tree.Name]:
            log.debug('%s %s' % (first, last))
            if first == 1:
                gls_list.append([(tree.Name, first)] + last)
            else:
                gls_list.append(last)

        log.debug("Number of inferred scenarios: %s" % len(d[tree.Name]))
        log.debug("Number of decisions: %s" % search_space)

        # the tracer stores all scores
        tracer = []

        for i, line in enumerate(gls_list):

            # calculate gains and losses
            gains = sum([1 for x in line if x[1] == 1])
            losses = sum([1 for x in line if x[1] == 0])

            # calculate the score
            if mode == 'w':
                score = r[0] * gains + r[1] * losses
            else:
                score = gains + losses

            # append it to the tracer
            tracer.append(score)

        # get the minimum score
        minScore = min(tracer)

        if mode == 'w':
            # return the minimal indices, sort them according to the number of
            # gains inferred, thereby pushing gains to the root, similar to
            # Mirkin's (2003) suggestion
            best_gls = [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore]
            best_gls = sorted(
                best_gls,
                key=lambda x: sum([i[1] for i in x]),
                reverse=push_gains
            )
            return best_gls[0]

            # push gains down to the root as suggested by Mirkin 2003
        minimal_gains = [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore]

        # make sure to check the model with minimal amount of gains
        minGains = len(self.taxa)
        for i, line in enumerate(minimal_gains):
            gains = sum([1 for x in line if x[1] == 1])
            if gains <= minGains:
                minGains = gains
        minimal_gains = [line for line in minimal_gains if sum(
            [1 for x in line if x[1] == 1]
        ) == minGains]

        best_scenario = 0
        old_length_of_tips = len(self.taxa) + 1

        for i, line in enumerate(minimal_gains):

            # calculate number of tips for the gains of a given scenario
            new_length_of_tips = 0
            for taxon, state in line:
                if state == 1:
                    new_length_of_tips += len(
                        self.tree.getNodeMatchingName(taxon).getTipNames()
                    )
            if new_length_of_tips < old_length_of_tips:
                old_length_of_tips = new_length_of_tips
                best_scenario = i

        return minimal_gains[best_scenario]

    def get_GLS(
        self,
        mode='weighted',
        ratio=(1, 1),
        restriction=3,
        output_gml=False,
        output_plot=False,
        tar=False,
        **keywords
    ):
        """
        Create gain-loss-scenarios for all non-singleton paps in the data.

        Parameters
        ----------
        mode : string (default="weighted")
            Select between "weighted", "restriction" and "topdown". The three
            modes refer to the following frameworks:

            * "weighted" refers to the weighted parsimony framework described in
              :evobib:`List2014b` and :evobib:`List2014a`.
              Weights are
              specified with help of a ratio for the scoring of gain and loss
              events. The ratio can be defined with help of the *ratio*
              keyword.
            * "restrictino" refers to a simple method in which only a
              specific amount of gain events
              is allowed. The maximally allowed number of gain events can be
              defined with help of the *restriction* keyword.
            * "topdown" refers to the top-down method outlined in
              :evobib:`Dagan2007` and first applied to linguistic data in
              :evobib:`Nelson-Sathi2011`. This method also defines a maximal
              number of gain events, but in contrast to the "restriction"
              approach, it starts from the top of the tree and stops if the
              maximal number of restrictions has been reached. The maximally
              allowed number of gain events can, again, be specified with help
              of the *restriction* keyword.
        ratio : tuple (default=(1,1))
            If "weighted" mode is selected, define the ratio between the
            weights for gains and losses.
        restriction : int (default=3)
            If "restriction" is selected as mode, define the maximal number of
            gains.
        output_gml : bool (default=False)
            If set to c{True}, the decisions for each GLS are stored in a
            separate file in GML-format.
        tar : bool (default=False)
            If set to c{True}, the GML-files will be added to a compressed tar-file.
        gpl : int (default=1)
            Specifies the maximal number of gains per lineage. This parameter
            specifies how cases should be handled in which a character is first
            gained, then lost, and then gained again. By setting this parameter
            to 1 (the default setting), such cases are prohibited, since only
            one gain per lineage is allowed.
        missing_data : int (default=0)
            Currently, we offer two ways to handle missing data. The first case
            just treats missing data in the same way in which the absence of a
            character is handled and can be evoked by setting this parameter to
            0. The second case will treat missing data as either absent or
            present characters, based on how well each option coincides with
            the overall evolutionary scenario. This behaviour can be evoked by
            setting this parameter to -1.
        push_gains: bool (default=True)
            In bottom-up calculations, there will often be multiple scenarios
            upon which only one is selected by the method. In order to define
            consistent criteria for scenario selection, we follow
            :evobib:`Mirkin2003` in allowing to force the algorithm to prefer
            those scenarios in which gains are pushed to the leaves. This
            behaviour is handle by this parameter. Setting it to *True* will
            force the algorithm to push gain events to the leaves of the tree.
            Setting it to *False* will force it to prefer those scenarios where
            the gains are closer to the root.


        """
        if mode not in ['weighted', 'w', 'r', 'restriction', 't', 'topdown']:
            raise ValueError("[!] The mode {0} is not available".format(mode))

        # define alias for mode
        if mode in ['w', 'weighted']:
            mode = 'weighted'
        elif mode in ['r', 'restriction']:
            mode = 'restriction'
        else:
            mode = 'topdown'

        # create a named string for the mode
        if mode == 'weighted':
            glm = 'w-{0[0]}-{0[1]}'.format(ratio)
        elif mode == 'restriction':
            glm = 'r-{0}'.format(restriction)
        elif mode == 'topdown':
            glm = 't-{0}'.format(restriction)

        # set defaults
        defaults = {
            "force": False,
            "gpl": 1,
            "push_gains": True,
            "missing_data": 0
        }
        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        # check for previous analyses
        if glm in self.gls and not keywords['force']:
            log.info(
                "Gain-loss scenario {0} has already been calculated.  For recalculation, "
                "set 'force' to 'True'.".format(glm))
            return

        # create statistics for this run
        self.stats[glm] = {}

        # store the statistics
        self.stats[glm]['mode'] = mode
        self.stats[glm]['dataset'] = self.dataset

        # attribute stores all gls for each cog
        self.gls[glm] = {}

        # make a temporary hash in order to decrease the number of calls to the
        # algorithm
        cogDict = {}

        skip, nonskip = 0, 0
        with util.ProgressBar('GAIN-LOSS-MAPPING ({0})'.format(glm), len(self.cogs)):
            for cog in self.cogs:
                # check whether cog has already been calculated
                cogTuple = tuple(self.paps[cog])
                if cogTuple in cogDict:
                    skip += 1
                    log.debug(
                        "Skipping already calculated pattern for COG {0}...".format(cog))
                    self.gls[glm][cog] = cogDict[cogTuple]
                else:
                    nonskip += 1
                    log.debug("Calculating GLS for COG {0}...".format(cog))

                    # check for singletons
                    if sum([x for x in self.paps[cog] if x == 1]) == 1:
                        gls = [(self.taxa[self.paps[cog].index(1)], 1)]
                    else:
                        if mode == 'weighted':
                            gls = get_gls(
                                self.paps[cog],
                                self.taxa,
                                self.tree,
                                gpl=keywords['gpl'],
                                weights=ratio,
                                push_gains=keywords['push_gains'],
                                missing_data=keywords['missing_data']
                            )

                        if mode == 'restriction':
                            gls = self._get_GLS(
                                self.paps[cog],
                                r=restriction,
                                mode='r',
                                gpl=keywords['gpl'],
                                push_gains=keywords['push_gains'],
                                missing_data=keywords['missing_data']
                            )

                        if mode == 'topdown':
                            gls = self._get_GLS_top_down(
                                self.paps[cog],
                                mode=restriction,
                                missing_data=keywords['missing_data']
                            )
                    noo = sum([t[1] for t in gls])

                    self.gls[glm][cog] = (gls, noo)

                    # append new results to cogDict
                    cogDict[cogTuple] = (gls, noo)

        # append scenario to gls
        log.info("Successfully calculated Gain-Loss-Scenarios.")

        # write the results to file
        # if output of gls is chosen, load the gml-graph
        if output_gml:
            # store the graph
            for cog in self.cogs:
                gls = self.gls[glm][cog][0]
                g = gls2gml(
                    gls,
                    self.tgraph,
                    self.tree,
                    filename=self._output_path('gml', '{0}-{1}'.format(self.dataset, glm),
                                               cog))

                # if plot of gml is chose
                if output_plot:
                    nodes = []

                    for n, d in g.nodes(data=True):
                        x = d['graphics']['x']
                        y = d['graphics']['y']
                        f = d['graphics']['fill']
                        o = d['origin']
                        l = d['label']

                        nodes.append((x, y, f, o, l))

                    edges = []
                    for a, b, d in g.edges(data=True):
                        xA = g.node[a]['graphics']['x']
                        xB = g.node[b]['graphics']['x']
                        yA = g.node[a]['graphics']['y']
                        yB = g.node[b]['graphics']['y']

                        edges += [(xA, xB, yA, yB)]

                    # mpl.rc('text',usetex=keywords['usetex'])
                    fig = plt.figure()
                    fig.add_subplot(111)
                    plt.axes(frameon=False)
                    plt.xticks([])
                    plt.yticks([])

                    plt.axis('equal')

                    for xA, xB, yA, yB in edges:
                        plt.plot(
                            [xA, xB],
                            [yA, yB],
                            '-',
                            color='black',
                            linewidth=5
                        )
                        plt.plot(
                            [xA, xB],
                            [yA, yB],
                            '-',
                            color='0.2',
                            linewidth=4
                        )
                    for x, y, f, o, l in nodes:
                        if f == '#000000':
                            c = '#ffffff'
                        else:
                            c = '#000000'
                        if o == 1:
                            size = 20
                        else:
                            size = 10
                        if l.startswith('edge') or l.startswith('root'):
                            plt.plot(x, y, 'o', markersize=size, color=f)
                        else:
                            plt.text(
                                x,
                                y,
                                l,
                                horizontalalignment='center',
                                verticalalignment='center',
                                size=8, fontweight='bold', color=c, backgroundcolor=f)

                    # plt.subplots_adjust(left=0.02,right=0.98,top=0.98,bottom=0.02)
                    plt.savefig(
                        self._output_path(
                            'gml', '{0}-{1}-figures'.format(self.dataset, glm), cog,
                            '.png'))
                    plt.clf()

            # if tar is chosen, put it into a tarfile
            if tar:
                os.system(
                    'cd {0}_phybo/gml/ ; tar -pczf {0}-{1}.tar.gz {0}-{1}; cd ..; cd ..'.format(
                        self.dataset,
                        glm
                    )
                )
                os.system('rm {0}_phybo/gml/{0}-{1}/*.gml'.format(self.dataset, glm))
                os.system('rmdir {0}_phybo/gml/{0}-{1}'.format(self.dataset, glm))

        # store some statistics as attributes
        self.stats[glm]['ano'] = sum(
            [v[1] for v in self.gls[glm].values()]
        ) / len(self.gls[glm])
        self.stats[glm]['mno'] = max([v[1] for v in self.gls[glm].values()])
        self.stats[glm]['ratio'] = ratio
        self.stats[glm]['restriction'] = restriction

        # store statistics and gain-loss-scenarios in textfiles
        log.info("Writing GLS data to file... ")

        lines = ['PAP\tGainLossScenario\tNumberOfOrigins']
        for cog in sorted(self.gls[glm]):
            gls, noo = self.gls[glm][cog]
            lines.append(
                "{0}\t".format(cog) + ','.join(
                    ["{0}:{1}".format(a, b) for a, b in gls]
                ) + '\t' + text_type(noo)
            )
        self._write_file(
            os.path.join('gls', '{0}-{1}.gls'.format(self.dataset, glm)), lines)

        # print out average number of origins
        log.info("Average Number of Origins: {0:.2f}".format(self.stats[glm]['ano']))

        # write statistics to stats file
        lines = [
            'Number of PAPs (total): {0}'.format(len(self.paps)),
            'Number of PAPs (non-singletons): {0}'.format(len(self.gls[glm])),
            'Number of Singletons: {0}'.format(len(self.singletons)),
            'Average Number of Origins: {0:.2f}'.format(self.stats[glm]['ano']),
            'Maximum Number of Origins: {0}'.format(self.stats[glm]['mno']),
            'Mode: {0}'.format(mode),
        ]
        if mode == 'weighted':
            lines.append('Ratio: {0[0]} / {0[1]}'.format(ratio))
        elif mode == 'restriction':
            lines.append('Restriction: {0}'.format(restriction))
        self._write_file(
            os.path.join('stats', '{0}-{1}'.format(self.dataset, glm)), lines)
        return

    def get_CVSD(
        self,
    ):
        """
        Calculate the Contemporary Vocabulary Size Distribution (CVSD).

        """
        # -># define taxa and concept as attribute for convenience
        # ->taxa = self.taxa
        # ->concepts = self.concept #XXX do we need this? XXX

        # -># calculate vocabulary size
        # ->forms = []
        # ->meanings = []
        # ->for taxon in taxa:
        # ->    f = [x for x in set(
        # ->        self.get_list(col=taxon,entry=self._pap_string,flat=True)
        # ->        ) if x in self.cogs
        # ->        ]
        # ->    m = set([x.split(':')[1] for x in f])
        # ->    forms += [len(f)]
        # ->    meanings += [len(m)]
        # ->
        # -># store the stuff as an attribute
        # ->self.dists['contemporary'] = [x for x,y in zip(forms,meanings)] # XXX
        dists = []
        for t in self.taxa:
            paps = sorted(set([p for p in self.get_list(
                taxa=t,
                entry=self._pap_string,
                flat=True
            ) if p not in self.singletons]))
            forms = len(paps)
            dists += [forms]  # / concepts]
        self.dists['contemporary'] = dists
        log.info("Calculated the distributions for contemporary taxa.")
        return

    def get_AVSD(
        self,
        glm,
        **keywords
    ):
        """
        Function retrieves all pap s for ancestor languages in a given tree.
        """
        # get keywords and defaults
        defaults = {
            'proto': False,
            'force': False,
        }
        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        # check for already calculated glm
        # check for previous analyses
        if glm in self.dists and not keywords['force'] and glm != 'mixed':
            log.info(
                "Gain-loss scenario {0} has already been calculated. For recalculation, "
                "set 'force' to 'True'.".format(glm))
            log.info("For recalculation, set 'force' to True.")
            return

        # get acs with help of utils
        acs, dst = get_acs(self, glm, **keywords)

        # append stuff to dist
        self.dists[glm] = dst

        # append stuff to acs
        self.acs[glm] = {}
        if keywords['proto']:
            paps = [self[k, self._pap_string] for k in self]
            protos = [self[k, keywords['proto']] for k in self]
            p2p = dict(zip(paps, protos))
            pap2protos = lambda x: p2p[x]
        else:
            pap2protos = lambda x: x  # dict(zip(paps,paps))

        for k, v in acs.items():
            self.acs[glm][k] = [(p, self.pap2con[p], pap2protos(p)) for p in v]

        # -># define concepts for convenience
        # ->concepts = self.concepts # XXX do we need this? XXX
        # ->
        # -># get all internal nodes, i.e. the nontips and also the root
        # ->nodes = ['root'] + sorted(
        # ->        [node.Name for node in self.tree.nontips()],
        # ->        key=lambda x: len(self.tree.getNodeMatchingName(x).tips()),
        # ->        reverse = True
        # ->        )

        # -># retrieve scenarios
        # ->tmp = sorted([(a,b,c) for a,(b,c) in self.gls[glm].items()])
        # ->cog_list = [t[0] for t in tmp]
        # ->gls_list = [t[1] for t in tmp]
        # ->noo_list = [t[2] for t in tmp]

        # -># create a list that stores the paps
        # ->paps = [[0 for i in range(len(nodes))] for j in range(len(cog_list))]

        # -># iterate and assign values
        # ->for i,cog in enumerate(cog_list):
        # ->
        # ->    # sort the respective gls
        # ->    gls = sorted(
        # ->            gls_list[i],
        # ->            key = lambda x: len(self.tree.getNodeMatchingName(x[0]).tips()),
        # ->            reverse = True
        # ->            )

        # ->    # retrieve the state of the root
        # ->    if gls[0][1] == 1 and gls[0][0] == 'root':
        # ->        state = 1
        # ->    else:
        # ->        state = 0

        # ->    # assign the state of the root to all nodes
        # ->    paps[i] = [state for node in nodes]

        # ->    # iterate over the gls and assign the respective values to all
        # ->    # children
        # ->    # XXX note that here we assume that missing data is coded as
        # ->    # 0, so this should probably be adapted XXX
        # ->    for name,event in gls:
        # ->        if event == 1:
        # ->            this_state = 1
        # ->        else:
        # ->            this_state = 0

        # ->        # get the subtree nodes
        # ->        sub_tree_nodes = [node.Name for node in
        # ->                self.tree.getNodeMatchingName(name).nontips()]

        # ->        # assign this state to all subtree nodes
        # ->        for node in sub_tree_nodes:
        # ->            paps[i][nodes.index(node)] = this_state

        # -># get number of forms and number of meanings
        # -># extract cogs instead of numbers, XXX this can actually be done in the
        # -># step before, it's just for testing at the moment
        # ->for i,cog in enumerate(cog_list):
        # ->    for j,t in enumerate(paps[i]):
        # ->        if t == 1:
        # ->            paps[i][j] = cog
        # ->        else:
        # ->            pass
        # ->
        # -># get forms and meanings
        # ->forms = []
        # ->meanings = []
        # ->for i in range(len(paps[0])):
        # ->    f = set([x[i] for x in paps if x[i] != 0])
        # ->    m = set([x[i].split(':')[1] for x in paps if x[i] != 0])
        # ->    forms += [len(f)]
        # ->    meanings += [len(m)]

        # -># store the number of forms as an attribute
        # ->self.dists[glm] = [x for x,y in zip(forms,meanings)] # XXX

        # -># store results of the analyses, that is, all paps for each ancestral
        # -># node
        # ->cogs = [k[self.header['pap']] for k in self._data.values()]

        # -># search for proto as keyword
        # ->if keywords['proto']:
        # ->    protos = [k[self.header[keywords['proto']]] for k in
        # ->            self._data.values()]
        # ->    cogs2proto = dict(zip(cogs,protos))
        # ->else:
        # ->    cogs2proto = dict(zip(cogs,cogs))

        # -># store data in acs attribute (ancestral cognate states)
        # ->self.acs[glm] = {}
        # ->for i,n in enumerate(nodes):
        # ->    for j,p in enumerate(paps):
        # ->        c = paps[j][i]
        # ->        if c != 0:
        # ->            m = self.pap2con[c]
        # ->            p = cogs2proto[c]

        # ->            if n != 'root':
        # ->                node = self.tree.getNodeMatchingName(n)
        # ->                node = n #''.join(
        # ->                        #[x for x in str(node) if x not in '";()'+"'"]
        # ->                        #)#.replace('(','').replace(')','').replace(',','-')
        # ->            else:
        # ->                node = n
        # ->
        # ->            try:
        # ->                self.acs[glm][node] += [(c,m,p)]
        # ->            except:
        # ->                self.acs[glm][node] = [(c,m,p)]

        log.info("Calculated the distributions for ancestral taxa.")
        return

    def plot_ACS(
        self,
        glm,
        **keywords
    ):
        """
        Plot a tree in which the node size correlates with the size of the ancestral node.
        """
        defaults = dict(
            scaler=0.1,
            degree=180,
            change=lambda x: 2.5 * x,
            figsize=(10, 5),
            colormap=mpl.cm.jet,
            colors=True
        )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        # check for the model
        if glm not in self.acs:
            self.get_ACS(glm, **keywords)

        # create a dictionary for all nodes
        node_dict = {}

        # iterate over contemporary taxa first
        vsizes = []
        tmp = {}
        for taxon in self.taxa:
            # get all cognates that are not singletongs
            cogs = sorted(set([x for x in self.get_list(
                col=taxon, flat=True, entry='pap') if x in self.singletons]))

            # count the number of paps
            node_dict[taxon] = dict(nodesize=len(cogs) * keywords['scaler'])

            vsizes += [len(cogs)]
            tmp[taxon] = len(cogs)

        # iterate over internal nodes now
        for a, b in [(x, y) for x, y in self.tree.getNodesDict().items() if
                     x not in self.taxa]:

            if a != 'root':
                node = a
                # node = str(b).replace(')','').replace('(','').replace(',','-')
            else:
                node = 'root'

            node_dict[a] = dict(nodesize=len(self.acs[glm][node]) * keywords['scaler'])
            tmp[a] = len(self.acs[glm][node])
            vsizes += [len(self.acs[glm][node])]

        # define a color-function
        if keywords['colors']:
            vsizes = sorted(set(vsizes))
            cfunc = np.array(np.linspace(10, 256, 245), dtype='int')
            for node in node_dict:
                node_dict[node]['nodecolor'] = mpl.colors.rgb2hex(
                    keywords['colormap'](
                        cfunc[int(tmp[node] * 244 / max(vsizes))]
                    )
                )

        # add the stuff to keywords
        keywords['node_dict'] = node_dict

        # check for filename in keywords
        if 'filename' not in keywords:
            keywords['filename'] = self._output_path(glm + '_acs')

        # plot the tree
        plot_tree(self.tree, no_labels=True, **keywords)

    def get_IVSD(
        self,
        output_gml=False,
        output_plot=False,
        tar=True,
        leading_model=False,
        mixed_threshold=0.0,
        evaluation='mwu',
        **keywords
    ):
        """
        Calculate VSD on the basis of each item.

        """
        kw = dict(
            fileformat='png',
            homoplasy=0.05,
        )
        kw.update(keywords)

        # assign concept dict
        mixed_concepts = {}

        # define concepts and taxa for convenience
        concepts = self.concepts
        taxa = self.taxa

        # get all internal nodes, i.e. the nontips and also the root
        nodes = ['root'] + sorted(
            [node.Name for node in self.tree.nontips()],
            key=lambda x: len(self.tree.getNodeMatchingName(x).tips()),
            reverse=True
        )

        # make dictionary that stores the best models for each cognate set
        best_models = {}

        # make array for all nodes and a dict for the scenarios
        all_avsd = [0 for node in nodes]
        scenarios = {}

        # iterate over concepts
        for concept in concepts:

            # get paps
            tmp = self.get_dict(row=concept, entry=self._pap_string)

            # add to list if value is missing
            for taxon in taxa:
                if taxon not in tmp:
                    tmp[taxon] = []

            # calculate distribution for contemporary taxa
            cvsd = [len([i for i in tmp[j] if i in self.cogs]) for j in taxa]

            # calculate ancestral dists, get all paps first
            pap_set = [i for i in set(
                self.get_list(
                    row=concept,
                    entry=self._pap_string,
                    flat=True
                )
            ) if i not in self.singletons]

            # get the models
            if leading_model:
                models = [leading_model] + sorted(
                    [k for k in self.gls.keys() if k != leading_model]
                )
            else:
                models = sorted(list(self.gls.keys()))

            models = [m for m in models if
                      m != 'mixed' and self._pvalues[m] >= mixed_threshold]

            # get the scenarios
            avsd_list = []
            for idx, glm in enumerate(models):

                tmp_list = []
                queue = ['root']
                while queue:

                    # get the parent
                    parent = queue.pop(0)

                    # get paps of parent
                    parent_paps = [p[0] for p in self.acs[glm][parent]]

                    # count number of paps
                    forms = [f for f in pap_set if
                             f in parent_paps]  # self.acs[glm][parent]]
                    tmp_list += [len(forms)]

                    children = self.tree.getNodeMatchingName(parent).Children
                    for child in children:
                        if child not in self.taxa:
                            queue += [child.Name]

                avsd_list += [[int(form - form * kw['homoplasy'] + 0.5) for
                               form in tmp_list]]

            # calculate best distribution, we can use averages for this
            # purpose, since it seems that the kruskalwallis test or
            # mannwhitneyu does not really apply to this kind of data with lots
            # of small numbers XXX
            zp_vsd = []
            cvsd_set = set(cvsd)
            for avsd in avsd_list:
                if len(cvsd_set) == 1 and set(avsd):
                    zp_vsd.append((0, 0.0))
                else:
                    if evaluation in ['mwu', 'mannwhitneyu']:
                        vsd = sp.stats.mstats.kruskalwallis(  # mannwhitneyu(
                            cvsd,
                            avsd,
                            # use_continuity=False
                        )
                        zp_vsd.append((vsd[0], vsd[1]))
                    elif evaluation in ['average']:
                        # check for best median and best average
                        ave_cvsd = sum(cvsd) / len(cvsd)
                        ave_avsd = sum(avsd) / len(avsd)

                        score = abs(ave_cvsd - ave_avsd)
                        zp_vsd.append((1, score))

            # extract p-values
            p_vsd = [p for z, p in zp_vsd]
            if evaluation in ['mwu', 'mannwhitneyu']:
                maxP = max(p_vsd)
            elif evaluation in ['average']:
                maxP = min(p_vsd)

            # check for threshold
            # if leading_model:
            #    if True: #maxP >= mixed_threshold:
            #        maxIdx = p_vsd.index(maxP)
            #        best_model = models[maxIdx]
            #    else:
            #        maxIdx = 0
            #        best_model = leading_model
            #        maxP = p_vsd[0]
            # else:
            maxIdx = p_vsd.index(maxP)
            best_model = models[maxIdx]

            for p in pap_set:
                gls, noo = self.gls[best_model][p]
                best_models[p] = (best_model, noo, maxP)
                scenarios[p] = (gls, noo)

            # add sum to general model XXX start here XXX
            all_avsd = [a + b for a, b in zip(avsd_list[maxIdx], all_avsd)]

            # add to concepts
            mixed_concepts[concept] = {}
            mixed_concepts[concept]['mixed'] = maxP

            for i, m in enumerate(models):
                mixed_concepts[concept][m] = p_vsd[i]

        self.best_models = best_models
        # print(sum([n for m,n,o in best_models.values()]) / len(best_models))

        # append to distributions

        # self.dists['mixed'] = all_avsd

        # append to available models
        self.gls['mixed'] = scenarios
        self.get_AVSD('mixed', **kw)

        # write the results to file
        # if output of gls is chosen, load the gml-graph
        if output_gml:
            # store the graph
            for cog in self.cogs:
                gls = self.gls["mixed"][cog][0]
                g = gls2gml(
                    gls,
                    self.tgraph,
                    self.tree,
                    filename=self._output_path(
                        'gml',
                        '{0}-{1}'.format(
                            self.dataset,
                            "mixed"
                        ),
                        cog
                    ),
                )

                # if plot of gml is chose
                if output_plot:
                    nodes = []

                    for n, d in g.nodes(data=True):
                        x = d['graphics']['x']
                        y = d['graphics']['y']
                        f = d['graphics']['fill']
                        o = d['origin']
                        l = d['label']

                        nodes.append((x, y, f, o, l))

                    edges = []
                    for a, b, d in g.edges(data=True):
                        xA = g.node[a]['graphics']['x']
                        xB = g.node[b]['graphics']['x']
                        yA = g.node[a]['graphics']['y']
                        yB = g.node[b]['graphics']['y']

                        edges += [(xA, xB, yA, yB)]

                    # mpl.rc('text',usetex=keywords['usetex'])
                    fig = plt.figure()
                    fig.add_subplot(111)
                    plt.axes(frameon=False)
                    plt.xticks([])
                    plt.yticks([])

                    plt.axis('equal')

                    for xA, xB, yA, yB in edges:
                        plt.plot(
                            [xA, xB],
                            [yA, yB],
                            '-',
                            color='black',
                            linewidth=5
                        )
                        plt.plot(
                            [xA, xB],
                            [yA, yB],
                            '-',
                            color='0.2',
                            linewidth=4
                        )
                    for x, y, f, o, l in nodes:
                        if f == '#000000':
                            c = '#ffffff'
                        else:
                            c = '#000000'
                        if o == 1:
                            size = 20
                        else:
                            size = 10
                        if l.startswith('edge') or l.startswith('root'):
                            plt.plot(x, y, 'o', markersize=size, color=f)
                        else:
                            plt.text(
                                x,
                                y,
                                l,
                                horizontalalignment='center',
                                verticalalignment='center',
                                size=8, fontweight='bold', color=c, backgroundcolor=f)

                    # plt.subplots_adjust(left=0.02,right=0.98,top=0.98,bottom=0.02)
                    plt.savefig(
                        self._output_path(
                            'gml',
                            '{0}-{1}-figures'.format(self.dataset, 'mixed'),
                            cog + '.' + kw['fileformat']))
                    plt.clf()

            # if tar is chosen, put it into a tarfile
            if tar:
                # FIXME: the code below is not portable to windows systems!
                os.system(
                    'cd {0}_phybo/gml/ ; tar -pczf {0}-{1}.tar.gz {0}-{1}; cd ..; cd ..'.format(
                        self.dataset,
                        "mixed"
                    )
                )
                os.system('rm {0}_phybo/gml/{0}-{1}/*.gml'.format(self.dataset, "mixed"))
                os.system('rmdir {0}_phybo/gml/{0}-{1}'.format(self.dataset, "mixed"))

        # store some statistics as attributes
        self.stats['mixed'] = {}
        self.stats['mode'] = 'mixed'
        self.stats['dataset'] = self.dataset
        self.stats['mixed']['ano'] = sum(
            [v[1] for v in self.gls['mixed'].values()]
        ) / len(self.gls['mixed'])
        self.stats['mixed']['mno'] = max([v[1] for v in self.gls['mixed'].values()])

        self.stats['mixed_concepts'] = mixed_concepts

        # store statistics and gain-loss-scenarios in textfiles
        log.info("Writing GLS data to file... ")

        # write gls-data to folder
        with util.TextFile(
            self._output_path('gls', '{0}-{1}.gls'.format(self.dataset, "mixed")),
            log=False
        ) as f:
            f.write('PAP\tGainLossScenario\tNumberOfOrigins\n')
            for cog in sorted(self.gls["mixed"]):
                gls, noo = self.gls["mixed"][cog]
                f.write(
                    "{0}\t".format(cog) + ','.join(
                        ["{0}:{1}".format(a, b) for a, b in gls]
                    ) + '\t' + str(noo) + '\n'
                )
        return

    def get_ACS(
        self,
        glm,
        **keywords
    ):
        """
        Compute the ancestral character states (ACS) for all internal nodes.

        """
        defaults = dict(
            proto=False,
            force=False,
            filename=self._output_path('acs-' + glm),
            fileformat='csv'
        )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        if glm not in self.acs:
            self.get_AVSD(glm, **keywords)
        elif keywords['force']:
            self.get_AVSD(glm, **keywords)

        with util.TextFile(keywords['filename'] + '.' + keywords['fileformat']) as f:
            for key in sorted(self.acs[glm].keys(), key=lambda x: len(x)):
                for c, m, p in sorted(self.acs[glm][key], key=lambda x: x[1]):
                    f.write('{0}\t{1}\t{2}\t{3}\n'.format(key, c, m, p))

    def get_MLN(
        self,
        glm,
        threshold=1,
        method='mr'
    ):
        """
        Compute an Minimal Lateral Network for a given model.

        Parameters
        ----------
        glm : str
            The dictionary key for the gain-loss-model.
        threshold : int (default=1)
            The threshold used to exclude edges.
        method  : str (default='mr')
            Select the method for MLN calculation. Choose between:
            * "mr": majority-rule, multiple links are resolved by selecting
              those which occur most frequently
            * "td": tree-distance, multiple links are resolved by selecting
              those which are closest on the tree
            * "bc": betweenness-centrality, multiple links are resolved by
              selecting those which have the highest betweenness centrality

        """

        # if not colormap:
        #    colormap = mpl.cm.jet

        # create the primary graph
        gPrm = nx.Graph()

        # make alias for tree and taxa for convenience
        taxa = self.taxa
        tree = self.tree

        # get the topological graph
        gTpl = self.tgraph

        # make alias for the current gls for convenience
        scenarios = self.gls[glm]

        # create dictionary for inferred lateral events
        ile = {}

        # create mst graph
        gMST = nx.Graph()

        # create out graph
        gOut = nx.Graph()

        # load data for nodes into new graph
        for node, data in gTpl.nodes(data=True):
            if data['label'] in taxa:
                data['graphics']['fill'] = '#ff0000'
                data['graphics']['type'] = 'rectangle'
                data['graphics']['w'] = 80.0
                data['graphics']['h'] = 20.0
            else:
                data['graphics']['type'] = 'ellipse'
                data['graphics']['w'] = 30.0
                data['graphics']['h'] = 30.0
                data['graphics']['fill'] = '#ff0000'
            gPrm.add_node(data['label'], **data)

        # load edge data into new graph
        for nodeA, nodeB, data in gTpl.edges(data=True):
            if 'graphics' not in data:
                data['graphics'] = {}
            data['graphics']['width'] = 10.0
            data['graphics']['fill'] = '#000000'
            data['label'] = 'vertical'

            gPrm.add_edge(
                gTpl.node[nodeA]['label'],
                gTpl.node[nodeB]['label'],
                **data
            )

        # start to assign the edge weights
        for cog, (gls, noo) in scenarios.items():

            # get the origins
            oris = [x[0] for x in gls if x[1] == 1]

            # connect origins by edges
            for i, oriA in enumerate(oris):
                for j, oriB in enumerate(oris):
                    if i < j:
                        try:
                            gPrm.edge[oriA][oriB]['weight'] += 1
                        except:
                            gPrm.add_edge(
                                oriA,
                                oriB,
                                weight=1
                            )

        log.info("Calculated primary graph.")
        log.info("Inferring lateral edges...")

        # create MST graph
        gMST = nx.Graph()

        with util.ProgressBar('MLN-REONSTRUCTION', len(scenarios)) as progress:
            for cog, (gls, noo) in scenarios.items():
                progress.update()
                ile[cog] = []

                # get the origins
                oris = [x[0] for x in gls if x[1] == 1]

                # create a graph of weights
                gWeights = nx.Graph()

                # calculate majority-rule edges
                if method in ['majority_rule', 'mr']:
                    # iterate over nodes
                    for i, nodeA in enumerate(oris):
                        for j, nodeB in enumerate(oris):
                            if i < j:
                                w = gPrm.edge[nodeA][nodeB]['weight']
                                gWeights.add_edge(
                                    nodeA,
                                    nodeB,
                                    weight=w
                                )
                elif method in ['tree_distance', 'td']:
                    for i, nodeA in enumerate(oris):
                        for j, nodeB in enumerate(oris):
                            if i < j:
                                try:
                                    w = len(
                                        self.tree.getConnectingEdges(
                                            nodeA,
                                            nodeB
                                        )
                                    )
                                except ValueError:
                                    if 'root' in (nodeA, nodeB):
                                        w = len(
                                            self.tree.getConnectingEdges(
                                                nodeB,
                                                nodeA
                                            )
                                        )
                                    else:
                                        wA = len(
                                            self.tree.getConnectingEdges(
                                                'root',
                                                nodeA
                                            )
                                        )
                                        wB = len(
                                            self.tree.getConnectingEdges(
                                                'root',
                                                nodeB
                                            )
                                        )
                                        w = wA + wB

                                gWeights.add_edge(
                                    nodeA,
                                    nodeB,
                                    weight=w
                                )
                elif method in ['betweenness_centrality', 'bc']:
                    bc = nx.edge_betweenness_centrality(
                        gPrm, normalized=True,
                        weight='weight'
                    )
                    for i, nodeA in enumerate(oris):
                        for j, nodeB in enumerate(oris):
                            if i < j:
                                try:
                                    w = bc[nodeA, nodeB]
                                except KeyError:
                                    w = bc[nodeB, nodeA]
                                # be careful with zero division
                                # if w == 0:
                                #    w = 0.1

                                gWeights.add_edge(
                                    nodeA,
                                    nodeB,
                                    weight=int(100 * (1 - w))  # int(1000 / w)
                                )
                elif method in ['central_node', 'cn']:

                    # get the weighted degrees for the primary graph
                    degrees = gPrm.degree(weight='weight')

                    # get the maximum degree
                    max_deg = sorted(
                        degrees,
                        key=lambda x: degrees[x],
                        reverse=True
                    )[0]

                    # add all nodes as simple
                    for i, nodeA in enumerate(oris):
                        for j, nodeB in enumerate(oris):
                            if i < j:
                                if max_deg in [nodeA, nodeB]:
                                    w = 0
                                else:
                                    w = 10

                                gWeights.add_edge(
                                    nodeA,
                                    nodeB,
                                    weight=w
                                )

                # if the graph is not empty
                if gWeights:

                    # check for identical weights and change them according to
                    # tree-distance
                    tmp_weights = {}
                    for a, b, d in gWeights.edges(data=True):
                        try:
                            tmp_weights[int(d['weight'])] += [(a, b)]
                        except:
                            tmp_weights[int(d['weight'])] = [(a, b)]

                    if method in ['mr', 'majority_rule']:
                        # check for identical weights and calculate the tree distance
                        for w in tmp_weights:
                            elist = tmp_weights[w]

                            # check whether there are more identical weights
                            if len(elist) > 1:

                                # if so, order all stuff according to branch
                                # length, we need try-except statement for
                                # branchdistances here, since cogent does not
                                # calculate distances to the root and back
                                branches = []
                                for a, b in elist:
                                    try:
                                        branch_distance = len(
                                            self.tree.getConnectingEdges(a, b))
                                        branches += [(a, b, branch_distance)]
                                    except:
                                        if 'root' in (a, b):
                                            branch_distance = len(
                                                self.tree.getConnectingEdges(b, a)
                                            )
                                            branches += [(a, b, branch_distance)]
                                        else:
                                            bdA = len(
                                                self.tree.getConnectingEdges('root', a)
                                            )
                                            bdB = len(
                                                self.tree.getConnectingEdges('root', b)
                                            )
                                            branches += [(a, b, bdA + bdB)]

                                # now change the weights according to the order
                                scaler = 1 / len(branches)
                                minus = 1 - scaler
                                branches = sorted(branches,
                                                  key=lambda x: (x[2], x[1], x[0]),
                                                  reverse=True)
                                for a, b, d in branches:
                                    gWeights.edge[a][b]['weight'] += minus
                                    minus -= scaler

                        # change maximum weights to distance weights
                        for a, b, d in sorted(gWeights.edges(data=True),
                                              key=lambda x: x[2]['weight']):
                            w = d['weight']
                            gWeights.edge[a][b]['weight'] = int(1000 / w) ** 2

                    # calculate the MST
                    mst = nx.minimum_spanning_tree(gWeights, weight='weight')

                    # assign the MST-weights to gMST
                    for nodeA, nodeB in mst.edges():
                        try:
                            gMST.edge[nodeA][nodeB]['weight'] += 1
                            gMST.edge[nodeA][nodeB]['cogs'] += [cog]
                        except:
                            gMST.add_edge(
                                nodeA,
                                nodeB,
                                weight=1,
                                cogs=[cog]
                            )
                        ile[cog] += [(nodeA, nodeB)]

        # load data for nodes into new graph
        for node, data in gTpl.nodes(data=True):
            if data['label'] in taxa:
                data['graphics']['fill'] = '#ff0000'
                data['graphics']['type'] = 'rectangle'
                data['graphics']['w'] = 80.0
                data['graphics']['h'] = 20.0
            else:
                data['graphics']['type'] = 'ellipse'
                data['graphics']['w'] = 30.0
                data['graphics']['h'] = 30.0
                data['graphics']['fill'] = '#ff0000'

            gOut.add_node(data['label'], **data)

        # load edge data into new graph
        for nodeA, nodeB, data in gTpl.edges(data=True):
            data['graphics']['width'] = 10.0
            data['graphics']['fill'] = '#000000'
            data['label'] = 'vertical'
            try:
                del data['graphics']['Line']
            except:
                pass
            gOut.add_edge(
                gTpl.node[nodeA]['label'],
                gTpl.node[nodeB]['label'],
                **data
            )

        # assign new edge weights
        for nodeA, nodeB, data in gMST.edges(data=True):
            w = data['weight']

            # get the color for the weight
            # color = mpl.colors.rgb2hex(colormap(cfunc[weights.index(w)]))

            data['graphics'] = {}
            # data['graphics']['fill'] = color
            # data['graphics']['width'] = w * scale
            data['cogs'] = ','.join([text_type(i) for i in data['cogs']])
            data['label'] = 'horizontal'

            # check for threshold
            if w >= threshold:
                try:
                    gOut.edge[nodeA][nodeB]
                except:
                    # add the data to the out-graph
                    gOut.add_edge(
                        nodeA,
                        nodeB,
                        **data
                    )
        # transfer node data

        log.info("Writing graph to file...")
        self._write_file('mln-' + glm + '.gml', nx.generate_gml(gOut))

        # write the inferred borrowing events (ILS, inferred lateral event)
        # between all taxa to file
        log.info("Writing Inferred Lateral Events to file...")

        with util.TextFile(self._output_path('ile-' + glm + '.csv')) as f:
            for cog, events in ile.items():
                if events:
                    f.write(
                        text_type(cog) + '\t' + ','.join(
                            ['{0}:{1}'.format(a, b) for a, b in events]
                        ) + '\n')

        # create file name for node labels (cytoscape output)
        lines = ["node.label (class=java.lang.String)"]
        for taxon in taxa:
            lines.append('{0} = {1}'.format(taxon, taxon))
        self._write_file('node.label.NA', lines)

        # add gOut to graphattributes
        self.graph[glm] = gOut

        # write stats to file
        # get the degree
        nodes = tree.getNodeNames()

        dgr, wdgr = [], []
        for taxon in nodes:
            horizontals = [g for g in gOut[taxon] if 'weight' in gOut[taxon][g]]

            dgr.append(len(horizontals))
            wdgr.append(sum([gOut[taxon][g]['weight'] for g in horizontals]))

        sorted_nodes = sorted(
            zip(nodes, dgr, wdgr),
            key=lambda x: x[1],
            reverse=True
        )
        lines = []
        for n, d, w in sorted_nodes:
            lines.append(
                '{0}\t{1}\t{2}\t{3}'.format(
                    n, text_type(tree.getNodeMatchingName(n)), d, w))
        self._write_file('taxa-' + glm + '.stats', lines)

        log.info("Wrote node degree distributions to file.")

        # write edge distributions
        edges = [g for g in gOut.edges(data=True) if 'weight' in g[2]]

        lines = []
        for nA, nB, d in sorted(
            edges,
            key=lambda x: x[2]['weight'],
            reverse=True
        ):
            lines.append(
                '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(
                    nA,
                    nB,
                    d['weight'],
                    d['cogs'],
                    tree.getNodeMatchingName(nA),
                    tree.getNodeMatchingName(nB)
                )
            )
        self._write_file('edge-' + glm + '.stats', lines)
        log.info("Wrote edge-weight distributions to file.")

        for taxon in self.taxa:
            with util.TextFile(
                self._output_path(os.path.join('taxa-' + glm, taxon + '.csv')), log=False
            ) as f:
                keys = [n for n in gOut[taxon] if gOut[taxon][n]['label'] == 'horizontal']
                for key in sorted(keys, key=lambda x: gOut[taxon][x]['weight']):
                    for cog in sorted(gOut[taxon][key]['cogs'].split(',')):
                        tmp = [x for x in self.etd[cog] if x != 0]
                        idx = [x[0] for x in tmp][0]
                        concept = self[idx, 'concept']

                        proto = cog

                        # get the index of the current entry in its dictionary
                        # representation
                        idx = self.get_dict(col=taxon, entry='pap')[concept]
                        idx = idx.index(cog)

                        # get its real index
                        idx = self.get_dict(col=taxon)[concept][idx]

                        # include entries specified in keywords XXX modify later
                        # for customization
                        for entry in ['ipa', 'proto']:
                            if entry in self.header:
                                proto += '\t' + self[idx, entry]

                        f.write('{0}\t{1}\t{2}\n'.format(key, proto, concept))
        log.info("Wrote list of edges per taxa to file.")
        return

    def get_PDC(
        self,
        glm,
        **keywords
    ):
        """
        Calculate Patchily Distributed Cognates.
        """
        defaults = dict(
            aligned_output=True,
        )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        patchy = {}
        paps = []

        for key, (gls, noo) in self.gls[glm].items():

            # get the origins
            oris = sorted(
                [x[0] for x in gls if x[1] == 1],
                key=lambda x: len(
                    self.tree.getNodeMatchingName(x).getTipNames()
                )
            )

            # get the tip-taxa for each origin
            tips = []

            losses = [a for a, b in zip(self.taxa, self.paps[key]) if b == 0]

            tipsofar = []
            for i, ori in enumerate(oris):
                new_tips = [
                    i + 1,
                    [t for t in
                     self.tree.getNodeMatchingName(
                         ori
                     ).getTipNames() if t not in losses and t not in tipsofar]
                ]
                tipsofar += new_tips[1]

                tips += [tuple(new_tips)]

            # now, all set of origins with their tips are there, we store them
            # in the patchy dictionary, where each taxon is assigned the
            # numerical value of the given patchy dist
            patchy[key] = {}
            if len(tips) > 1:
                for i, tip in tips:
                    for taxon in tip:
                        patchy[key][taxon] = i
            else:
                for i, tip in tips:
                    for taxon in tip:
                        patchy[key][taxon] = 0

            paps.append((key, noo))

        log.info("[i] Retrieved patchy distributions.")

        # get the index for the paps in the wordlist
        papIdx = self.header['pap']
        taxIdx = self._colIdx

        # create a dictionary as updater for the wordlist
        updater = {}
        for key in self:

            # get the taxon first
            taxon = self[key][taxIdx]

            # get the pap
            pap = self[key][papIdx]

            try:
                updater[key] = '{0}:{1}'.format(pap, patchy[pap][taxon])
            except KeyError:
                updater[key] = '{0}:{1}'.format(pap, 0)

        # update the wordlist
        self.add_entries(
            'patchy',
            updater,
            lambda x: x,
            override=True
        )

        # write data to file
        # self.output('csv',filename=self.dataset+'_phybo/wl-'+glm)
        # XXX change later

        log.info("[i] Updated the wordlist.")

        # write ranking of concepts to file
        with util.TextFile(self._output_path('paps-' + glm + '.stats')) as f:
            if 'proto' in self.entries:
                f.write('COGID\tGLID\tCONCEPT\tORIGINS\tREFLEXES\tORIG/REFL\tPROTO\n')
            else:
                f.write('COGID\tGLID\tCONCEPT\tORIGINS\tREFLEXES\tORIG/REFL\n')
            concepts = {}
            for a, b in sorted(paps, key=lambda x: x[1], reverse=True):
                a1, a2 = a.split(':')
                a3 = self._id2gl[int(a2)]

                # check for number of occurrences
                l = [k for k in self.etd[a] if k != 0]

                # append three vals: number of origins, number of words, and the
                # number of origins per number of words
                try:
                    concepts[a3] += [(b, len(l), b / len(l))]
                except:
                    concepts[a3] = [(b, len(l), b / len(l))]

                # check for proto
                if 'proto' in self.entries:
                    proto = self[[k[0] for k in l][0], 'proto']
                    f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2f}\t{6}\n'.format(
                        a1, a2, a3, b, len(l), b / float(len(l)), proto))
                else:
                    f.write(
                        '{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2f}\n'.format(a1, a2, a3, b, len(l),
                                                                    float(b) / len(l)))
        log.info("Wrote stats on paps to file.")

        # write stats on concepts
        with util.TextFile(self._output_path('concepts-' + glm + '.stats')) as f:
            f.write('CONCEPT\tORIGINS\tREFLEXES\tORI/REF\n')

            cstats = {}
            for key in concepts:
                # get origins per concept
                oriperc = sum([c[0] for c in concepts[key]]) / len(concepts[key])

                # get patchy cognate sets per number of words
                patchyperw = sum([c[2] for c in concepts[key]]) / len(concepts[key])

                # get the number of words per concept
                numperc = sum([c[1] for c in concepts[key]]) / len(concepts[key])

                cstats[key] = (oriperc, numperc, patchyperw)
                # concepts[key] = sum(concepts[key])/len(concepts[key])

            for a, b in sorted(cstats.items(), key=lambda x: x[1][2], reverse=True):
                f.write('{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n'.format(a, b[0], b[1], b[2]))

            # write average to file
            f.write('TOTAL\t{0:.2f}\t{1:.2f}\t{2:.2f}\n'.format(
                sum([cstats[c][0] for c in cstats]) / len(cstats),
                sum([cstats[c][1] for c in cstats]) / len(cstats),
                sum([cstats[c][2] for c in cstats]) / len(cstats)
            ))
        log.info("Wrote stats on concepts to file.")

        # write alternative stats on concepts including information of
        # singletons (excluding them may bias the results)
        with util.TextFile(self._output_path('cognates-' + glm + '.stats')) as f:
            f.write('CONCEPT\tCOGNATES\tPATCHIES\tREFLEXES\tPCR\n')

            concepts = {}
            for pap in self.etd:
                gloss = self.pap2con[pap]
                idxs = [idx[0] for idx in self.etd[pap] if idx != 0]
                patchies = [self[idx, 'patchy'] for idx in idxs]
                cogs = [self[idx, self._pap_string] for idx in idxs]
                reflexes = len(patchies)
                patchies = len(set(patchies))
                cogs = len(set(cogs))
                try:
                    concepts[gloss] += [(cogs, patchies, reflexes)]
                except:
                    concepts[gloss] = [(cogs, patchies, reflexes)]

            for key, value in concepts.items():
                concepts[key] = (
                    sum([v[0] for v in value]),
                    sum([v[1] for v in value]),
                    sum([v[2] for v in value])
                )

            for k, (c, p, r) in sorted(concepts.items(), key=lambda x: x[1][2],
                                       reverse=True):
                f.write(
                    '{0}\t{1}\t{2}\t{3}\t{4:.2f}\n'.format(k, c, p, r, (p - c + 1) / r))
            # write mean
            mc = sum([x[0] for x in concepts.values()]) / len(concepts)
            mp = sum([x[1] for x in concepts.values()]) / len(concepts)
            mr = sum([x[2] for x in concepts.values()]) / len(concepts)

            f.write('{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\t{4:.2f}\n'.format(
                "MEAN", mc, mp, mr, (mp - mc + 1) / mr))

        log.info("Wrote stats on concepts to file.")

        # store params in attribute stats
        self.stats["CONCEPTS"] = dict(
            cognates=sum([x[0] for x in concepts.values()]) / self.height,
            patchies=sum([x[1] for x in concepts.values()]) / self.height,
            reflexes=sum([x[2] for x in concepts.values()]) / self.height,
            origins=sum([x[0] for x in cstats.values()]) / self.height,
            # reflexes = sum([x[1] for x in cstats.values()]) / self.height,
            patchy_per_reflex=sum([x[2] for x in cstats.values()]) / self.height
        )

        # write results to alm-file
        # get all patchy cognates
        tmp = {}
        for key in self:
            patchy = self[key, 'patchy']
            if not patchy.endswith('0'):

                concept = self[key, 'concept']
                taxon = self[key, 'doculect']
                pap = self[key, 'pap']

                # XXX change this later for more flexibility XXX

                if 'ipa' in self.header:
                    word = self[key, 'ipa']
                else:
                    word = self[key, 'counterpart']
                if not word:
                    raise NameError(
                        "[ERROR] Neither 'ipa' nor 'counterpart' is defined."
                    )

                if concept not in tmp:
                    tmp[concept] = {}
                if pap not in tmp[concept]:
                    tmp[concept][pap] = {}

                try:
                    tmp[concept][pap][patchy] += [(taxon, word)]
                except:
                    tmp[concept][pap][patchy] = [(taxon, word)]

        if keywords["aligned_output"]:
            # write stuff to alm-file
            with util.TextFile(
                self._output_path(self.dataset + '-' + glm + '.alm.patchy'), log=False
            ) as f:
                for concept in sorted(tmp.keys()):
                    f.write('# Basic Concept: "{0}"\n\n'.format(concept))

                    for pap in sorted(tmp[concept].keys()):
                        f.write('## Cognate-Set: "{0}"\n'.format(pap))

                        words = []
                        langs = []
                        patchies = []

                        for patchy in sorted(tmp[concept][pap].keys()):
                            # get words and languages
                            words += [t[1].replace("ˈ", '') for t in
                                      tmp[concept][pap][patchy]]
                            langs += [t[0] for t in tmp[concept][pap][patchy]]

                            patchies += [patchy[-1] for i in
                                         range(len(tmp[concept][pap][patchy]))]

                        msa = Multiple(words)
                        # XXX add for different alignment algorithm later XXX
                        msa.prog_align()
                        alms = msa.alm_matrix

                        # get formatter for languages
                        formatter = max([len(lang) for lang in langs])

                        for i, word in enumerate(words):
                            string = '{0:' + text_type(
                                formatter) + '}\t{1}\t|\t{2}\t|\t[{3}]\n'
                            f.write(
                                string.format(langs[i], patchies[i], '\t'.join(alms[i]),
                                              word))
                        f.write('\n')
                    f.write('\n')

    def get_edge(
        self,
        glm,
        nodeA,
        nodeB,
        entries='',
        msn=False
    ):
        """
        Return the edge data for a given gain-loss model.
        """
        # define a warning message
        warning = "No edge between {0} and {1} could be found".format(
            nodeA,
            nodeB
        )
        # check for entryB
        if isinstance(entries, text_type):
            entries = entries.split(',')

        # get the graph locally for convenience
        if not msn:
            graph = self.graph[glm]
        else:
            graph = self.geograph[glm]

        # get the edge
        try:
            edge = graph.edge[nodeA][nodeB]
        except:
            log.info(warning)
            return

        # check the edge
        if not msn:
            if edge['label'] == 'horizontal':
                cogs = edge['cogs'].split(',')
            else:
                log.info(warning)
                return
        else:
            cogs = edge['cogs'].split(',')

        # define list for output
        outA = {}
        outB = {}

        # check whether nodes are in list or not
        if nodeA in self.taxa:
            nodesA = [nodeA]
        else:
            nodesA = self.tree.getNodeMatchingName(nodeA).getTipNames()
        if nodeB in self.taxa:
            nodesB = [nodeB]
        else:
            nodesB = self.tree.getNodeMatchingName(nodeB).getTipNames()

        # assemble the data
        outA = {}
        for node in nodesA:
            tmp = dict(
                zip(
                    self.get_list(
                        col=node,
                        flat=True,
                        entry='pap'
                    ),
                    self.get_list(
                        col=node,
                        flat=True,
                    )
                )
            )
            for cog in cogs:
                vals = [node]
                for entry in entries:
                    try:
                        vals += [self[tmp[cog], entry]]
                    except:
                        pass
                if len(vals) > 1:
                    try:
                        outA[cog] += [tuple(vals)]
                    except:
                        outA[cog] = [tuple(vals)]

        # assemble the data
        outB = {}
        for node in nodesB:
            tmp = dict(
                zip(
                    self.get_list(
                        col=node,
                        flat=True,
                        entry='pap'
                    ),
                    self.get_list(
                        col=node,
                        flat=True,
                    )
                )
            )
            for cog in cogs:
                vals = [node]
                for entry in entries:
                    try:
                        vals += [self[tmp[cog], entry]]
                    except:
                        pass
                if len(vals) > 1:
                    try:
                        outB[cog] += [tuple(vals)]
                    except:
                        outB[cog] = [tuple(vals)]

        # assemble the output
        output = []
        for cog in cogs:
            try:
                output += [
                    (
                        self.pap2con[cog], outA[cog], outB[cog]
                    )
                ]
            except:
                self.log.error(
                    "Error encountered in cognate {0}.".format(self.pap2con[cog]))
        return output

    def analyze(
        self,
        runs="default",
        mixed=False,
        output_gml=False,
        tar=False,
        full_analysis=True,
        plot_dists=False,
        output_plot=False,
        plot_mln=False,
        plot_msn=False,
        **keywords
    ):
        """
        Carry out a full analysis using various parameters.

        Parameters
        ----------
        runs : {str list} (default="default")
            Define a couple of different models to be analyzed. Select between:

            * 'default': weighted analysis, using parsimony and weights for
              gains and losses
            * 'topdown': use the traditional approach by
              :evobib:`Nelson-Sathi2011`
            * 'restriction': use the restriction approach

            You can also define your own mix of models.

        usetex : bool (default=True)
            Specify whether you want to use LaTeX to render plots.
        mixed : bool (default=False)
            If set to c{True}, calculate a mixed model by selecting the best
            model for each item separately.
        output_gml : bool (default=False)
            Set to c{True} in order to output every gain-loss-scenario in
            GML-format.
        full_analysis : bool (default=True)
            Specifies whether a full analysis is carried out or not.
        plot_mln : bool (default=True)
            Select or unselect output plot for the MLN.
        plot_msn : bool (default=False)
            Select or unselect output plot for the MSN.

        """

        # set defaults
        defaults = {
            "colorbar": None,  # mpl.cm.jet,
            'threshold': 1,
            'fileformat': rcParams['phybo_fileformat'],
            'usetex': False,
            'only': [],
            'colormap': None,  # mpl.cm.jet
            'proto': False,
            'xticksize': 6,
            'method': 'mr',  # majority rule
            'gpl': 1,
            "push_gains": True,
            "missing_data": 0,
            "aligned_output": False,
            "homoplasy": 0.05,
            'evaluation': 'mwu'
        }

        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        # define a default set of runs
        if runs in ['default', 'weighted']:
            runs = [
                ('weighted', (3, 1)),
                ('weighted', (5, 2)),
                ('weighted', (2, 1)),
                ('weighted', (3, 2)),
                ('weighted', (1, 1)),
            ]

        elif runs in ['topdown', 'top-down']:
            runs = [('topdown', 2),
                    ('topdown', 3),
                    ('topdown', 4),
                    ('topdown', 5),
                    ('topdown', 6),
                    ('topdown', 7),
                    ('topdown', 8),
                    ('topdown', 9),
                    ('topdown', 10),
                    ]

        elif runs == 'restriction':

            runs = [('restriction', 2),
                    ('restriction', 3),
                    ('restriction', 4),
                    ('restriction', 5),
                    ('restriction', 6),
                    ]

        # carry out the various analyses
        for mode, params in runs:
            if mode == 'weighted':
                log.info(
                    "Analysing dataset with mode {0} ".format(mode) +
                    "and ratio {0[0]}:{0[1]}...".format(params)
                )

                self.get_GLS(
                    mode=mode,
                    ratio=params,
                    output_gml=output_gml,
                    tar=tar,
                    output_plot=output_plot,
                    gpl=keywords['gpl'],
                    push_gains=keywords['push_gains'],
                    missing_data=keywords["missing_data"],
                )
            elif mode == 'restriction':
                log.info(
                    "Analysing dataset with mode {0} ".format(mode) +
                    "and restriction {0}...".format(params)
                )

                self.get_GLS(
                    mode=mode,
                    restriction=params,
                    output_gml=output_gml,
                    tar=tar,
                    output_plot=output_plot,
                    gpl=keywords['gpl'],
                    push_gains=keywords['push_gains'],
                    missing_data=keywords["missing_data"]
                )
            elif mode == 'topdown':
                log.info(
                    "Analysing dataset with mode {0} ".format(mode) +
                    "and restriction {0}...".format(params)
                )
                self.get_GLS(
                    mode=mode,
                    restriction=params,
                    output_gml=output_gml,
                    tar=tar,
                    output_plot=output_plot,
                    missing_data=keywords["missing_data"]
                )

        # calculate the different distributions
        # start by calculating the contemporary distributions
        log.info("Calculating the Contemporary Vocabulary Distributions...")
        self.get_CVSD()

        # now calculate the rest of the distributions
        log.info("Calculating the Ancestral Vocabulary Distributions...")

        modes = list(self.gls.keys())
        with util.ProgressBar('ANCESTRAL VOCABULARY DISTRIBUTIONS',
                              len(modes)) as progress:
            for m in modes:
                progress.update()
                self.get_AVSD(m, **keywords)

        # compare the distributions using mannwhitneyu
        log.info("Comparing the distributions...")

        zp_vsd = []
        for m in modes:
            vsd = sp.stats.mstats.kruskalwallis(  # sp.stats.mannwhitneyu(
                self.dists['contemporary'],
                self.dists[m],
                # use_continuity = False
            )

            zp_vsd.append((vsd[0], vsd[1]))

        # determine the best model
        p_vsd = [p for z, p in zp_vsd]
        maxP = max(p_vsd)
        glm = modes[p_vsd.index(maxP)]

        # make a hash for all ps for all models
        self._pvalues = dict(zip(modes, p_vsd))

        # set the best model
        self.best_model = glm

        # calculate mixed model
        if mixed:
            log.info("Calculating the mixed model...")
            self.get_IVSD(
                output_plot=output_plot,
                output_gml=output_gml,
                tar=tar,
                leading_model=glm,
                **keywords
            )

            # set the mixed model as the best one
            self.best_model = 'mixed'

            if 'mixed' not in modes:
                modes += ['mixed']
                vsd = sp.stats.mstats.kruskalwallis(  # sp.stats.mannwhitneyu(
                    self.dists['contemporary'],
                    self.dists['mixed'],
                    # use_continuity = False
                )
                zp_vsd.append((vsd[0], vsd[1]))

        # write results to file
        log.info("Writing stats to file.")
        with util.TextFile(self._output_path(self.dataset + '.stats')) as f:
            f.write("Mode\tANO\tMNO\tVSD_z\tVSD_p\n")
            for i in range(len(zp_vsd)):
                f.write(
                    '{0}\t{1:.2f}\t{2}\t{3}\n'.format(
                        modes[i],
                        self.stats[modes[i]]['ano'],
                        self.stats[modes[i]]['mno'],
                        '{0[0]}\t{0[1]:.4f}'.format(zp_vsd[i])
                    )
                )

        # plot the stats if this is defined in the settings
        if plot_dists:

            log.info("Plotting distributions.")
            # specify latex
            mpl.rc('text', usetex=keywords['usetex'])

            # store distributions in lists
            dists_vsd = [self.dists[m] for m in modes]

            # store contemporary dists
            dist_vsd = self.dists['contemporary']

            # get the average number of origins
            ano = [self.stats[m]['ano'] for m in modes]

            # create a sorter for the distributions
            sorter = [s[0] for s in sorted(
                zip(range(len(modes)), ano),
                key=lambda x: x[1]
            )]

            # sort the stuff
            dists_vsd = [dists_vsd[i] for i in sorter]
            modes = [modes[i] for i in sorter]
            mode_strings = [m for m in modes]

            # sort the zp-values
            zp_vsd = [zp_vsd[i] for i in sorter]

            # format the zp-values
            if keywords['usetex']:

                p_vsd = []
                for i, (z, p) in enumerate(zp_vsd):
                    if p < 0.001:
                        p_vsd.append('p$<${0:.2f}'.format(p))
                    elif p >= 0.05:
                        p_vsd.append(r'\textbf{{p$=${0:.2f}}}'.format(p))

                        # adjust the modes
                        mode_strings[i] = r'\textbf{' + modes[i] + '}'
                    else:
                        p_vsd.append('p$=${0:.2f}'.format(p))

            else:
                p_vsd = []
                for z, p in zp_vsd:
                    if p < 0.001:
                        p_vsd.append('p<{0:.2f}'.format(p))
                    elif p >= 0.05:
                        p_vsd.append(r'p={0:.2f}'.format(p))
                    else:
                        p_vsd.append('p={0:.2f}'.format(p))

            # create the figure
            fig = plt.figure()

            # create the axis
            ax = fig.add_subplot(111)

            # add the boxplots
            b = ax.boxplot([dist_vsd] + dists_vsd)
            plt.setp(b['medians'], color='black')
            plt.setp(b['whiskers'], color='black')
            plt.setp(b['boxes'], color='black')

            # adjust the yticks
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(18)

            # add the xticks
            plt.xticks(
                range(1, len(modes) + 2),
                [''] + ['{0}\n{1}'.format(m, p) for m, p in zip(mode_strings, p_vsd)],
                size=keywords['xticksize'],
                rotation=45,
                # rotation_mode='anchor',
                ha='center',
                # va = 'center'
            )

            ax.yaxis.grid(
                True,
                linestyle='-',
                which='major',
                color='lightgrey',
                alpha=0.5,
                zorder=1
            )

            plt.subplots_adjust(bottom=0.2)

            # save the figure
            plt.savefig(self._output_path('vsd.' + keywords['fileformat']))
            plt.clf()

            log.info("Plotted the distributions.")

        # carry out further analyses if this is specified
        if full_analysis:

            self.get_MLN(
                self.best_model,
                threshold=keywords['threshold'],
                method=keywords['method']
            )

            # check whether plots are chosen
            if plot_mln:
                self.plot_MLN(
                    self.best_model,
                    filename=self._output_path('mln-' + glm),
                    threshold=keywords['threshold'],
                    fileformat=keywords['fileformat'],
                    usetex=keywords['usetex'],
                    colormap=keywords['colormap']
                )
            if plot_msn:
                self.plot_MSN(
                    self.best_model,
                    filename=self._output_path('msn-' + glm),
                    fileformat=keywords['fileformat'],
                    threshold=keywords['threshold'],
                    only=keywords['only'],
                    usetex=keywords['usetex'],
                    colormap=keywords['colormap']
                )

            self.get_PDC(
                self.best_model,
                **keywords
            )

    def plot_MLN(
        self,
        glm='',
        fileformat='pdf',
        threshold=1,
        usetex=False,
        taxon_labels='taxon_short_labels',
        alphat=False,
        alpha=0.75,
        **keywords
    ):
        """
        Plot the MLN with help of Matplotlib.

        glm : str (default='')
            Identifier for the gain-loss model that is plotted. Defaults to the
            model that had the best scores in terms of probability.
        filename : str (default='')
            If no filename is selected, the filename is identical with the
            dataset.
        fileformat : {'svg','png','jpg','pdf'} (default='pdf')
            Select the format of the output plot.
        threshold : int (default=1)
            Select the threshold for drawing lateral edges.
        usetex : bool (default=True)
            Specify whether you want to use LaTeX to render plots.
        colormap : {None matplotlib.cm}
            A :py:class:`matplotlib.colormap` instance. If set to c{None}, this
            defaults to :py:class:`~matplotlib.cm.jet`.
        taxon_labels : str (default='taxon.short_labels')
            Specify the taxon labels that should be included in the plot.

        """
        # check for correct glm
        if not glm and hasattr(self, 'best_model'):
            glm = self.best_model
        elif not glm:
            raise ValueError(
                "[i] You should select an appropriate model first."
            )

        # switch backend, depending on whether tex is used or not
        backend = mpl.get_backend()
        if usetex and backend != 'pgf':
            plt.switch_backend('pgf')
            mpl.rcParams['text.latex.unicode'] = True
        elif not usetex and backend != 'TkAgg':
            plt.switch_backend('TkAgg')

        defaults = dict(
            figsize="optimal",  # rcParams['phybo_figsize'],
            figure_width=10,
            figure_scale=1,
            colormap=mpl.cm.jet,
            filename=self.dataset,
            linescale=rcParams['phybo_linescale'],
            maxweight=rcParams['phybo_maxweight'],
            xlim=rcParams['phybo_xlim'],
            ylim=rcParams['phybo_ylim'],
            xlimr=rcParams['phybo_xlimr'],
            xliml=rcParams['phybo_xliml'],
            ylimt=rcParams['phybo_ylimt'],
            ylimb=rcParams['phybo_ylimb'],
            left=rcParams['phybo_left'],
            right=rcParams['phybo_right'],
            top=rcParams['phybo_top'],
            bottom=rcParams['phybo_bottom'],
            cbar_shrink=rcParams['phybo_cbar_shrink'],
            cbar_fraction=rcParams['phybo_cbar_fraction'],
            cbar_pad=rcParams['phybo_cbar_pad'],
            cbar_orientation=rcParams['phybo_cbar_orientation'],
            cbar_label=rcParams['phybo_cbar_label'],
            vedgestyle=rcParams['phybo_vedgestyle'],
            vedgecolor=rcParams['phybo_vedgecolor'],
            vedgelinewidth=rcParams['phybo_vedgelinewidth'],
            vedgeinnerline=rcParams['phybo_vedgeinnerline'],
            hedgescale=rcParams['phybo_hedgescale'],
            nodestyle=rcParams['phybo_nodestyle'],
            nodesize=rcParams['phybo_nodesize'],
            nodecolor=rcParams['phybo_nodecolor'],
            labels=rcParams['phybo_labels'],
            _prefix=rcParams['phybo_prefix'],
            _suffix=rcParams['phybo_suffix'],
            textsize=rcParams['phybo_textsize'],
            vsd_scale=rcParams['phybo_vsd_scale'],
            latex_preamble=rcParams['phybo_latex_preamble'],
            alpha_threshold=0.2
        )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        # get max and min values for coordinates
        xvals, yvals = [], []
        for n, d in self.graph[glm].nodes(data=True):
            xvals += [d['graphics']['x']]
            yvals += [d['graphics']['y']]

        minX, minY, maxX, maxY = min(xvals), min(yvals), max(xvals), max(yvals)
        if keywords['figsize'] == 'optimal':
            w = maxX + abs(minX)
            h = maxY + abs(minY)
            keywords['figsize'] = (
                keywords['figure_width'] + keywords['figure_scale'],
                h / (w / (keywords['figure_width'])
                     )
            )

        if keywords['latex_preamble']:
            mpl.rcParams['pgf.preamble'] = keywords['latex_preamble']

        colormap = keywords['colormap']
        filename = keywords['filename']

        # define labels
        labels = {}
        for taxon in self.taxa:
            if taxon not in keywords['labels']:
                labels[taxon] = taxon
            else:
                labels[taxon] = keywords['labels'][taxon]

        # create a dictionary for all nodes
        node_dict = {}

        # iterate over contemporary taxa first
        for taxon in self.taxa:
            # get all cognates that are not singletongs
            cogs = [x for x in
                    self.get_list(col=taxon, flat=True, entry='pap') if x in self.cogs]

            # count the number of paps
            node_dict[taxon] = len(cogs) * keywords['vsd_scale']

        # iterate over internal nodes now
        for a, b in [(x, y) for x, y in self.tree.getNodesDict().items() if
                     x not in self.taxa]:

            if a != 'root':
                node = a

            else:
                node = 'root'

            if node in self.acs[glm]:
                node_dict[a] = len(self.acs[glm][node]) * keywords['vsd_scale']
            else:
                node_dict[a] = 1 * keywords['vsd_scale']

        # get the graph
        graph = self.graph[glm]

        # store in internal and external nodes
        inodes = []
        enodes = []

        # get colormap for edgeweights
        edge_weights = []
        for nodeA, nodeB, data in graph.edges(data=True):
            if data['label'] == 'horizontal':
                edge_weights.append(data['weight'])

        # add max weight to edge_weights
        if keywords['maxweight']:
            edge_weights += range(keywords['maxweight'])

        # determine a colorfunction
        cfunc = np.array(np.linspace(10, 256, len(set(edge_weights))), dtype='int')

        # sort the weights
        weights = sorted(set(edge_weights))

        # get the scale for the weights (needed for the line-width)
        try:
            scale = 20.0 / max(edge_weights) * keywords['linescale']
        except ValueError:
            scale = 1.0

        # write colors and scale to graph
        for nA, nB, data in graph.edges(data=True):
            if data['label'] == 'horizontal':
                w = data['weight']
                data['graphics'] = {}
                data['graphics']['fill'] = mpl.colors.rgb2hex(
                    colormap(cfunc[weights.index(w)]))
                data['graphics']['width'] = scale * w

        # get the nodes
        for n, d in graph.nodes(data=True):
            g = d['graphics']
            x = g['x']
            y = g['y']
            h = g['h']
            w = g['w']
            s = g['s']

            # get the nodesize
            if keywords['nodestyle'] == 'vsd':
                try:
                    ns = node_dict[n]
                except:
                    ns = keywords['nodesize']
            else:
                ns = keywords['nodesize']

            if d['label'] not in self.taxa:
                inodes += [(x, y, ns)]
            else:
                if 'angle' in d['graphics']:
                    r = d['graphics']['angle']
                else:
                    r = 0

                enodes += [(x, y, d['label'], r, s, ns)]

        # store vertical and lateral edges
        vedges = []
        ledges = []
        weights = []

        # get the edges
        for a, b, d in graph.edges(data=True):

            xA = graph.node[a]['graphics']['x']
            yA = graph.node[a]['graphics']['y']
            xB = graph.node[b]['graphics']['x']
            yB = graph.node[b]['graphics']['y']

            if d['label'] == 'vertical':

                vedges += [(xA, xB, yA, yB)]
            else:
                g = d['graphics']
                f = g['fill']
                w = g['width']
                a = alpha
                if d['weight'] < threshold:
                    if alphat:
                        a = keywords['alpha_threshold']
                    else:
                        w = 0.0

                ledges += [(xA, xB, yA, yB, f, w, a)]

                weights.append(d['weight'])

        if not weights:
            weights = [0]

        # usetex
        mpl.rc('text', usetex=usetex)

        # create the figure
        fig = plt.figure(
            facecolor='white',
            figsize=keywords['figsize']
        )
        figsp = fig.add_subplot(111)

        # create the axis
        plt.axes(frameon=False)
        plt.xticks([0], [''])
        plt.yticks([0], [''])

        # set equal axis
        plt.axis('equal')

        # draw the horizontal edges
        for xA, xB, yA, yB, f, w, a in sorted(ledges, key=lambda x: x[-2]):
            plt.plot(
                [xA, xB],
                [yA, yB],
                '-',
                color=f,
                linewidth=float(w) / keywords['hedgescale'],
                alpha=a
            )

        # draw the vertical edges
        for xA, xB, yA, yB in vedges:
            plt.plot(
                [xA, xB],
                [yA, yB],
                '-',
                color=keywords['vedgecolor'],
                linewidth=keywords['vedgelinewidth'],
            )
            if keywords['vedgestyle'] == 'double':
                plt.plot(
                    [xA, xB],
                    [yA, yB],
                    '-',
                    color='1.0',
                    linewidth=keywords['vedgeinnerline'],
                )
        # store x,y values for ylim,xlim drawing
        xvals = []
        yvals = []

        # draw the nodes
        for x, y, s in inodes:
            xvals += [x]
            yvals += [y]

            plt.plot(
                x,
                y,
                'o',
                markersize=s,  # keywords['nodesize'],
                color=keywords['nodecolor'],
            )
            if keywords['nodestyle'] == 'double':
                plt.plot(
                    x,
                    y,
                    'o',
                    markersize=s,  # keywords['nodesize']-4,
                    color='white'
                )

        for x, y, t, r, ha, s in enodes:

            xvals += [x]
            yvals += [y]

            # plot the marker
            plt.plot(
                x,
                y,
                'o',
                markersize=s,  # keywords['nodesize'],
                color=keywords['nodecolor'],
                zorder=200
            )

            if keywords['nodestyle'] == 'double':
                plt.plot(
                    x,
                    y,
                    'o',
                    markersize=s,  # keywords['nodesize']-4,
                    color='white'
                )

            # this is a workaround to get the text away from the node
            if ha == 'left':
                text = keywords['_prefix'] + labels[t]
            else:
                text = labels[t] + keywords['_suffix']

            plt.text(
                x,
                y,
                text,
                size=keywords['textsize'],
                verticalalignment='center',
                horizontalalignment=ha,
                fontweight='bold',
                color='black',
                rotation=r,
                rotation_mode='anchor',
                zorder=1
            )

        # add a colorbar
        cax = figsp.imshow(
            [[1, 2], [1, 2]],
            cmap=colormap,
            visible=False
        )
        cbar = fig.colorbar(
            cax,
            ticks=[
                1,
                1.25,
                1.5,
                1.75,
                2
            ],
            orientation=keywords['cbar_orientation'],
            shrink=keywords['cbar_shrink'],
            fraction=keywords['cbar_fraction'],
            pad=keywords['cbar_pad']
        )

        # check for maxweights-keyword
        if keywords['maxweight']:
            weights += [keywords['maxweight']]

        cbar.set_clim(1.0)
        cbar.set_label('Inferred Links')
        cbar.ax.set_yticklabels(
            [
                text_type(min(weights)),
                '',
                text_type(int(max(weights) / 2)),
                '',
                text_type(max(weights))
            ]
        )

        if keywords['xliml'] and keywords['xlimr']:
            xliml = keywords['xliml']
            xlimr = keywords['xlimr']
        else:
            xliml, xlimr = 2 * [keywords['xlim']]

        if keywords['ylimb'] and keywords['ylimt']:
            ylimb = keywords['ylimb']
            ylimt = keywords['ylimt']
        else:
            ylimb, ylimt = 2 * [keywords['ylim']]

        plt.ylim(min(yvals) - ylimb, max(yvals) + ylimt)
        plt.xlim(min(xvals) - xliml, max(xvals) + xlimr)
        plt.subplots_adjust(
            left=keywords['left'],
            right=keywords['right'],
            top=keywords['top'],
            bottom=keywords['bottom']
        )
        # fig.axes.get_xaxis().set_visible(False)
        # fig.axes.get_yaxis().set_visible(False)
        plt.axis('off')

        # save the figure
        plt.savefig(filename + '.' + fileformat)
        plt.clf()
        log.file_written(filename + '.' + fileformat)
        return

    def plot_MLN_3d(
        self,
        glm='',
        filename='',
        fileformat='pdf',
        threshold=1,
        usetex=True,
        colormap=None,  # mpl.cm.jet,
        taxon_labels='taxon_short_labels',
        alphat=False,
        alpha=0.75,
        **keywords
    ):
        """
        Plot the MLN with help of Matplotlib in 3d.

        glm : str (default='')
            Identifier for the gain-loss model that is plotted. Defaults to the
            model that had the best scores in terms of probability.
        filename : str (default='')
            If no filename is selected, the filename is identical with the
            dataset.
        fileformat : {'svg','png','jpg','pdf'} (default='pdf')
            Select the format of the output plot.
        threshold : int (default=1)
            Select the threshold for drawing lateral edges.
        usetex : bool (default=True)
            Specify whether you want to use LaTeX to render plots.
        colormap : {None matplotlib.cm}
            A :py:class:`matplotlib.colormap` instance. If set to c{None}, this
            defaults to :py:class:`~matplotlib.cm.jet`.
        taxon_labels : str (default='taxon.short_labels')
            Specify the taxon labels that should be included in the plot.

        """
        # add azim and the other params
        kw = dict(
            azim=220,
            elev=22
        )
        kw.update(keywords)

        # check for correct glm
        if not glm and hasattr(self, 'best_model'):
            glm = self.best_model
        elif not glm:
            raise ValueError(
                "[i] You should select an appropriate model first."
            )

        # switch backend, depending on whether tex is used or not
        backend = mpl.get_backend()
        if usetex and backend != 'pgf':
            plt.switch_backend('pgf')
            mpl.rcParams['text.latex.unicode'] = True

        elif not usetex and backend != 'TkAgg':
            plt.switch_backend('TkAgg')

        # check for filename
        if not filename:
            filename = self.dataset

        # if not colormap
        if not colormap:
            colormap = mpl.cm.jet

        # set default, XXX change later
        if 'height' not in keywords:
            keywords['height'] = 7
        if 'width' not in keywords:
            keywords['width'] = 10

        # try to load the configuration file
        conf = self._config()

        # check for 'taxon.labels' in conf
        if taxon_labels in conf:  # XXX change later
            tfunc = lambda x: conf[taxon_labels][x]
        else:
            tfunc = lambda x: x

        # get the graph
        graph = self.graph[glm]

        # store in internal and external nodes
        inodes = []
        enodes = []

        # get colormap for edgeweights
        edge_weights = []
        for nodeA, nodeB, data in graph.edges(data=True):
            if data['label'] == 'horizontal':
                edge_weights.append(data['weight'])

        # determine a colorfunction
        cfunc = np.array(np.linspace(10, 256, len(set(edge_weights))), dtype='int')

        # sort the weights
        weights = sorted(set(edge_weights))

        # get the scale for the weights (needed for the line-width)
        scale = 20.0 / max(edge_weights)

        # write colors and scale to graph
        for nA, nB, data in graph.edges(data=True):
            if data['label'] == 'horizontal':
                w = data['weight']
                data['graphics'] = {}
                data['graphics']['fill'] = mpl.colors.rgb2hex(
                    colormap(cfunc[weights.index(w)]))
                data['graphics']['width'] = scale * w

        # get the nodes
        for n, d in graph.nodes(data=True):
            g = d['graphics']
            x = g['x']
            y = g['y']
            z = g['z']
            s = g['s']
            zorder = g['zorder']

            if d['label'] not in self.taxa:
                inodes += [(x, y, -z, zorder)]
            else:
                if 'angle' in d['graphics']:
                    r = d['graphics']['angle']
                else:
                    r = 0

                if usetex:
                    enodes += [(
                        x,
                        y,
                        -z,
                        r'\textbf{' + tfunc(d['label']).replace('_', r'\_') + r'}',
                        r,
                        s,
                        zorder
                    )]
                else:
                    enodes += [(x, y, -z, tfunc(d['label']), r, s, zorder)]

        # store vertical and lateral edges
        vedges = []
        ledges = []
        weights = []

        # get the edges
        for a, b, d in graph.edges(data=True):

            xA = graph.node[a]['graphics']['x']
            yA = graph.node[a]['graphics']['y']
            zA = graph.node[a]['graphics']['z']
            xB = graph.node[b]['graphics']['x']
            yB = graph.node[b]['graphics']['y']
            zB = graph.node[b]['graphics']['z']
            zoA = graph.node[a]['graphics']['zorder']
            zoB = graph.node[b]['graphics']['zorder']
            zorder = int((zoA + zoB) / 2)

            if d['label'] == 'vertical':
                vedges += [(xA, xB, yA, yB, -zA, -zB, zorder)]
            else:
                g = d['graphics']
                f = g['fill']
                w = g['width']
                a = alpha
                if d['weight'] < threshold:
                    if alphat:
                        a = 0.2
                    else:
                        w = 0.0

                ledges += [(xA, xB, yA, yB, -zA, -zB, zorder, f, w, a)]

                weights.append(d['weight'])

        # usetex
        mpl.rc('text', usetex=usetex)

        # create the figure
        fig = plt.figure(
            facecolor='white',
            figsize=(keywords['width'], keywords['height'])
        )
        figsp = fig.add_subplot(111, projection='3d')

        # draw the horizontal edges
        for xA, xB, yA, yB, zA, zB, zorder, f, w, a in sorted(ledges,
                                                              key=lambda x: x[-2]):
            figsp.plot(
                [xA, xB],
                [yA, yB],
                [zA, zB],
                color=f,
                linewidth=float(w) / 4,
                alpha=a,
                zorder=zorder  # 100 # * abs(xA-xB) + 100 * abs(yA-yB)
            )

        # draw the vertical edges
        for xA, xB, yA, yB, zA, zB, zorder in vedges:
            figsp.plot(
                [xA, xB],
                [yA, yB],
                [zA, zB],
                color='0.0',
                linewidth=3,
                alpha=0.75,
                zorder=zorder  # 100 * abs(xA-xB) + 100 * abs(yA-yB)
            )
            # figsp.plot(
            #        [xA,xB],
            #        [yA,yB],
            #        [zA,zB],
            #        color='1.0',
            #        linewidth=1,
            #        #zorder = 100 * abs(xA-xB) + 100 * abs(yA-yB)
            #        )
        # store x,y values for ylim,xlim drawing
        xvals = []
        yvals = []
        zvals = []

        # draw the nodes
        for x, y, z, zorder in inodes:
            xvals += [x]
            yvals += [y]
            zvals += [z]

            figsp.scatter(
                x,
                y,
                z,
                marker='o',
                s=20,
                c='black',
                zorder=zorder  # 100 * x + 100 * y
            )

        # draw the leaves
        # store x and y-maxima for ylim, xlim drawing
        for x, y, z, t, r, ha, zorder in enodes:
            xvals += [x]
            yvals += [y]
            zvals += [z]

            figsp.text(
                x,
                y,
                z,
                t,
                size='5',
                verticalalignment='center',
                horizontalalignment='center',
                bbox=dict(
                    facecolor='white',
                    boxstyle='square,pad=0.2',
                    ec="none",
                    # alpha = 0.25
                ),
                fontweight='bold',
                color='black',
                zorder=zorder + 200  # 120 # * x + 100 * y
            )

        figsp.view_init(azim=kw['azim'], elev=kw['elev'])
        figsp.set_ylim(min(yvals), max(yvals))
        figsp.set_xlim(min(xvals), max(xvals))
        figsp.set_zlim(min(zvals), max(zvals))

        figsp.set_axis_off()

        plt.savefig(filename + '.' + fileformat, bbbox_inches='tight')
        plt.clf()

    def get_MSN(
        self,
        glm='',
        fileformat='pdf',
        external_edges=False,
        deep_nodes=False,
        **keywords
    ):
        """
        Plot the Minimal Spatial Network.

        Parameters
        ----------
        glm : str (default='')
            A string that encodes which model should be plotted.
        filename : str
            The name of the file to which the plot shall be written.
        fileformat : str
            The output format of the plot.
        threshold : int (default=1)
            The threshold for the minimal amount of shared links that shall be
            plotted.
        usetex : bool (default=True)
            Specify whether LaTeX shall be used for the plot.

        """
        # check for correct glm
        if not glm and hasattr(self, 'best_model'):
            glm = self.best_model
        elif not glm:
            raise ValueError(
                "[i] You should select an appropriate model first."
            )

        # redefine taxa and tree for convenience
        taxa, tree = self.taxa, self.tree

        # get the graph
        graph = self.graph[glm]

        # XXX check for coordinates of the taxa, otherwise load them from file and
        # add them to the wordlist XXX add later, we first load it from file
        if 'coords' in self._meta:
            coords = self._meta['coords']
        else:
            coords = csv2dict(
                self.dataset,
                'coords',
                dtype=[str, float, float]
            )

        # calculate all resulting edges, using convex hull as
        # approximation
        geoGraph = nx.Graph()
        for node, data in graph.nodes(data=True):
            geoGraph.add_node(node, **data)

        edges = list(graph.edges(data=True))
        with util.ProgressBar('MINIMAL SPATIAL NETWORK', len(edges)) as progress:
            for nA, nB, d in edges:
                progress.update()

                # get the labels
                lA = graph.node[nA]['label']
                lB = graph.node[nB]['label']

                # first check, whether edge is horizontal
                if d['label'] == 'horizontal':

                    # if both labels occur in taxa, it is simple
                    if lA in taxa and lB in taxa:
                        try:
                            geoGraph.edge[lA][lB]['weight'] += d['weight']
                            geoGraph.edge[lA][lB]['cogs'] += ',' + d['cogs']
                        except:
                            geoGraph.add_edge(lA, lB, weight=d['weight'], cogs=d['cogs'])
                    elif not external_edges:
                        # if only one in taxa, we need the convex hull for that node
                        if lA in taxa or lB in taxa:

                            # check which node is in taxa
                            if lA in taxa:
                                this_label = lA
                                other_nodes = tree.getNodeMatchingName(lB).getTipNames()
                                other_label = lB
                            elif lB in taxa:
                                this_label = lB
                                other_nodes = tree.getNodeMatchingName(lA).getTipNames()
                                other_label = lA

                            # first, get all the cogs
                            cogs = d['cogs'].split(',')

                            # iterate over all cogs now
                            for cog in cogs:

                                # check whether the nodes have the respective cognate
                                # and take only those that have it
                                new_other_nodes = []
                                for other_node in other_nodes:
                                    paps = self.get_list(
                                        col=other_node,
                                        entry='pap',
                                        flat=True
                                    )
                                    if cog in paps and other_node != this_label:
                                        new_other_nodes += [other_node]

                                # get the convex points of others
                                these_coords = [
                                    (round(coords[t][0], 5), round(coords[t][1], 5)) for t
                                    in
                                    new_other_nodes]
                                hulls = getConvexHull(these_coords, polygon=False)

                                # get the hull with the minimal euclidean distance
                                distances = []
                                for hull in hulls:
                                    distances.append(np.linalg.norm(
                                        np.array(hull) - np.array(coords[this_label])))
                                this_hull = hulls[distances.index(min(distances))]
                                other_label = new_other_nodes[
                                    these_coords.index(
                                        (
                                            round(this_hull[0], 5),
                                            round(this_hull[1], 5)
                                        )
                                    )
                                ]

                                # append the edge to the graph
                                try:
                                    geoGraph.edge[this_label][other_label]['weight'] += 1
                                    geoGraph.edge[this_label][other_label][
                                        'cogs'] += ',' + cog
                                except:
                                    geoGraph.add_edge(this_label, other_label, weight=1,
                                                      cogs=cog)

                        elif deep_nodes:
                            # get the taxa of a and b
                            taxA = tree.getNodeMatchingName(lA).getTipNames()
                            taxB = tree.getNodeMatchingName(lB).getTipNames()

                            # get the cogs
                            cogs = d['cogs'].split(',')

                            # iterate over the cogs
                            for cog in cogs:
                                newtaxA = []
                                newtaxB = []

                                # get the lists
                                for t in taxA:
                                    paps = self.get_list(
                                        col=t,
                                        entry='pap',
                                        flat=False
                                    )
                                    if cog in paps:
                                        newtaxA += [t]
                                for t in taxB:
                                    paps = self.get_list(
                                        col=t,
                                        entry='pap',
                                        flat=False
                                    )
                                    if cog in paps:
                                        newtaxB += [t]

                                # get the convex points
                                coordsA = [
                                    (round(coords[t][0], 5), round(coords[t][1], 5)) for t
                                    in newtaxA]
                                coordsB = [
                                    (round(coords[t][0], 5), round(coords[t][1], 5)) for t
                                    in newtaxB]
                                hullsA = getConvexHull(coordsA, polygon=False)
                                hullsB = getConvexHull(coordsB, polygon=False)

                                # get the closest points
                                distances = []
                                hulls = []
                                for i, hullA in enumerate(hullsA):
                                    for j, hullB in enumerate(hullsB):
                                        distances.append(np.linalg.norm(
                                            np.array(hullA) - np.array(hullB)))
                                        hulls.append((hullA, hullB))
                                minHulls = hulls[distances.index(min(distances))]

                                labelA = newtaxA[coordsA.index(
                                    (round(minHulls[0][0], 5), round(minHulls[0][1], 5)))]
                                labelB = newtaxB[coordsB.index(
                                    (round(minHulls[1][0], 5), round(minHulls[1][1], 5)))]

                                # append the edge to the graph
                                try:
                                    geoGraph.edge[labelA][labelB]['weight'] += 1
                                    geoGraph.edge[labelA][labelB]['cogs'] += ',' + cog
                                except:
                                    geoGraph.add_edge(labelA, labelB, weight=1, cogs=cog)

        # write stats to file
        with util.TextFile(self._output_path('taxa-msn-' + glm + '.stats')) as f:
            # get the degree
            nodes = tree.getTipNames()

            dgr, wdgr = [], []
            for taxon in nodes:
                horizontals = [g for g in geoGraph[taxon] if
                               'weight' in geoGraph[taxon][g]]

                dgr.append(len(horizontals))
                wdgr.append(sum([geoGraph[taxon][g]['weight'] for g in horizontals]))

            sorted_nodes = sorted(
                zip(nodes, dgr, wdgr),
                key=lambda x: x[1],
                reverse=True
            )
            for n, d, w in sorted_nodes:
                f.write(
                    '{0}\t{1}\t{2}\t{3}\n'.format(
                        n, text_type(tree.getNodeMatchingName(n)), d, w))

        # write edge distributions
        with util.TextFile(self._output_path('edge-msn-' + glm + '.stats')) as f:
            edges = [g for g in geoGraph.edges(data=True) if 'weight' in g[2]]

            for nA, nB, d in sorted(
                edges,
                key=lambda x: x[2]['weight'],
                reverse=True
            ):
                f.write(
                    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                        nA,
                        nB,
                        d['weight'],
                        d['cogs'],
                        tree.getNodeMatchingName(nA),
                        tree.getNodeMatchingName(nB)
                    )
                )

        try:
            self.geograph[glm] = geoGraph
        except:
            self.geograph = {}
            self.geograph[glm] = geoGraph
        return

    def plot_MSN(
        self,
        glm='',
        fileformat='pdf',
        threshold=1,
        usetex=False,
        alphat=False,
        alpha=0.75,
        only=[],
        **keywords
    ):
        """
        Plot a minimal spatial network.
        """
        # set defaults
        defaults = dict(
            latex_preamble=[],
            figsize=(10, 10),
            colormap=mpl.cm.jet,
            filename=self.dataset,
            linescale=1.0,
            maxweight=False,
            xlim=5,
            ylim=5,
            xlimr=False,
            xliml=False,
            ylimt=False,
            ylimb=False,
            left=0.02,
            right=0.98,
            top=1.00,
            bottom=0.00,
            cbar_shrink=0.55,
            cbar_fraction=0.1,
            cbar_pad=0.1,
            cbar_orientation='vertical',
            cbar_label='Inferred Links',
            cbar_fontsize=10,
            resolution='l',
            table_text_color='black',
            water_color='0.2',
            lw=2,
            cmap_max=250,
            continent_color='0.9',
            projection='merc',
            legend_size=18,
            linewidth=4,
            min_lon=False,
            max_lon=False,
            min_lat=False,
            max_lat=False,
            table_column_width=[0.025, 0.1325],
            coastline_color="0.5",
            table_location=3,
            legend_location=[0.85, 0.02],
            table_cell_height=0.024,
            table_text_size=10,
            alpha=0.75,
            cmap_min=30,
            markersize=20
        )

        # load the rc-file XXX add internal loading later
        conf = self._config()
        if not conf:
            try:
                conf = self._meta['conf']
            except:
                raise ValueError('[!] Configuration is not specified!')

        # overwrite configuration from keywords
        for k in keywords:
            conf[k] = keywords[k]

        # overwrite keywords with defaults
        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        # check for only
        if not only:
            only = self.taxa

        # switch backend, depending on whether tex is used or not
        backend = mpl.get_backend()
        if usetex and backend != 'pgf':
            plt.switch_backend('pgf')
            mpl.rcParams['text.latex.unicode'] = True
        elif not usetex and backend != 'TkAgg':
            plt.switch_backend('TkAgg')

        # check for preamble settings
        if keywords['latex_preamble']:
            mpl.rcParams['pgf.preamble'] = keywords['latex_preamble']

        # usetex
        mpl.rc('text', usetex=usetex)

        # define stuff for convenience
        filename = keywords['filename']
        colormap = keywords['colormap']

        # check for groups, add functionality for groups in qlc-file later XXX
        if 'groups' in self._meta:
            groups = self._meta['groups']
        else:
            groups = dict([(k, v) for k, v in csv2list(self.dataset, 'groups')])

        # update configuration
        for k in keywords:
            if k not in conf:
                conf[k] = keywords[k]

        # set the graph variable
        geoGraph = self.geograph[glm]

        # get the weights for the lines
        weights = []
        for a, b, d in geoGraph.edges(data=True):
            weights += [d['weight']]
        max_weight = max(weights)
        sorted_weights = sorted(set(weights))

        # get a color-function
        color_dict = np.array(
            np.linspace(
                0,
                256,
                len(set(weights))
            ),
            dtype='int'
        )

        # get a line-function
        line_dict = np.linspace(
            0.5,
            conf['linewidth'],
            len(set(weights))
        )

        # XXX check for coordinates of the taxa, otherwise load them from file and
        # add them to the wordlist XXX add later, we first load it from file
        if 'coords' in self._meta:
            coords = self._meta['coords']
        else:
            coords = csv2dict(
                self.dataset,
                'coords',
                dtype=[str, float, float]
            )

        # determine the maxima of the coordinates
        latitudes = [i[0] for i in coords.values()]
        longitudes = [i[1] for i in coords.values()]

        min_lat, max_lat = min(latitudes), max(latitudes)
        min_lon, max_lon = min(longitudes), max(longitudes)

        # start to initialize the basemap
        fig = plt.figure(figsize=keywords['figsize'])
        figsp = fig.add_subplot(111)

        # instantiate the basemap
        m = bmp.Basemap(
            llcrnrlon=min_lon + conf['min_lon'],
            llcrnrlat=min_lat + conf['min_lat'],
            urcrnrlon=max_lon + conf['max_lon'],
            urcrnrlat=max_lat + conf['max_lat'],
            resolution=conf['resolution'],
            projection=conf['projection']
        )

        # draw first values
        m.drawmapboundary(fill_color=conf['water_color'])
        m.drawcoastlines(color=conf['continent_color'], linewidth=0.5)
        m.drawcountries(color=conf['coastline_color'], linewidth=0.5)
        m.fillcontinents(color=conf['continent_color'], lake_color=conf['water_color'])

        # plot the lines
        for a, b, d in sorted(geoGraph.edges(data=True), key=lambda x: x[2]['weight']):

            if a in coords and b in coords and a in only or b in only:
                w = d['weight']

                color = colormap(color_dict[sorted_weights.index(w)])
                linewidth = line_dict[sorted_weights.index(w)]

                if w < threshold:
                    if alphat:
                        alpha = 0.2
                    else:
                        linewidth = 0.0
                else:
                    alpha = conf['alpha']

                # retrieve the coords
                yA, xA = coords[a]
                yB, xB = coords[b]

                # get the points on the map
                xA, yA = m(xA, yA)
                xB, yB = m(xB, yB)

                # plot the points
                plt.plot(
                    [xA, xB],
                    [yA, yB],
                    '-',
                    color=color,
                    alpha=alpha,
                    linewidth=linewidth,
                    zorder=w + 50
                )

        # plot the points for the languages
        cell_text = []
        legend_check = []

        # check for taxon.labels in conf
        if 'taxon_labels' in conf:
            tfunc = lambda x: conf['taxon_labels'][x]
        else:
            tfunc = lambda x: x
        if 'groups_labels' in conf:
            gfunc = lambda x: conf['groups_labels'][x]
        else:
            gfunc = lambda x: x

        # check for defaults
        defaults = {
            "markersize": 10,
            "table_cell_height": 0.025,
        }
        for k in defaults:
            if k not in conf:
                conf[k] = defaults[k]

        for i, (taxon, (lng, lat)) in enumerate(
                sorted(coords.items(), key=lambda x: x[0])):

            # retrieve x and y from the map
            x, y = m(lat, lng)

            # get the color of the given taxon
            # taxon_color = colors[groups[taxon]]

            # get colors from conf
            this_group = groups[taxon]
            taxon_color = conf['groups_colors'][this_group]
            try:
                taxon_marker = conf['groups_markers'][this_group]
            except:
                taxon_marker = 'o'

            # check for legend

            if gfunc(groups[taxon]) in legend_check:
                # plot the marker
                plt.plot(
                    x,
                    y,
                    taxon_marker,
                    markersize=conf['markersize'],
                    color=taxon_color,
                    zorder=max_weight + 52,
                )
            else:
                # plot the marker
                plt.plot(
                    x,
                    y,
                    taxon_marker,
                    markersize=conf['markersize'],
                    color=taxon_color,
                    zorder=max_weight + 52,
                    label=gfunc(groups[taxon])
                )
                legend_check.append(gfunc(groups[taxon]))

            # add number to celltext
            if usetex:
                cell_text.append([str(i + 1), tfunc(taxon).replace('_', r'\_')])
            else:
                cell_text.append([str(i + 1), tfunc(taxon)])

            # plot the text
            # check for darkness of color
            if taxon_color in ['black', 'gray'] or taxon_color[:3] in ['0.3', '0.2',
                                                                       '0.1', '0.0']:
                text_color = 'white'
            else:
                text_color = 'black'

            plt.text(
                x,
                y,
                str(i + 1),
                size=str(int(conf['markersize'] / 2)),
                color=text_color,
                label=taxon,
                horizontalalignment='center',
                fontweight="bold",
                verticalalignment='center',
                zorder=max_weight + 55
            )

        # add a colorbar
        cax = figsp.imshow(
            [[1, 2], [1, 2]],
            visible=False,
            cmap=colormap
        )
        cbar = fig.colorbar(
            cax,
            ticks=[
                1,
                1.25,
                1.5,
                1.75,
                2
            ],
            orientation=keywords['cbar_orientation'],
            shrink=keywords['cbar_shrink'],
            fraction=keywords['cbar_fraction'],
            pad=keywords['cbar_pad'],
        )
        cbar.set_clim(1.0)
        cbar.set_label(keywords['cbar_label'], fontsize=keywords['cbar_fontsize'])
        cbar.ax.set_yticklabels(
            [
                str(min(weights)),
                '',
                str(int(max(weights) / 2)),
                '',
                str(max(weights))
            ],
            fontsize=keywords['cbar_fontsize']
        )

        # add the legend
        this_table = plt.table(
            cellText=cell_text,
            colWidths=conf['table_column_width'],
            loc=conf['table_location'],
        )
        this_table.auto_set_font_size(False)
        this_table.set_fontsize(conf['table_text_size'])

        # adjust the table
        for line in this_table._cells:
            this_table._cells[line]._text._horizontalalignment = 'left'
            this_table._cells[line]._text._fontproperties.set_weight('bold')
            this_table._cells[line]._text.set_color(conf['table_text_color'])
            this_table._cells[line].set_height(conf['table_cell_height'])
            # this_table._cells[line]._text._fontproperties.set_size(conf['table.text.size'])
            this_table._cells[line].set_linewidth(0.0)
            this_table._cells[line].set_color(conf['table_cell_color'])

        this_table.set_zorder(100)

        plt.legend(
            loc=conf['legend_location'],
            numpoints=1,
            prop={
                'size': conf['legend_size'],
                'weight': 'bold'
            }
        )

        plt.subplots_adjust(
            left=keywords['left'],
            right=keywords['right'],
            top=keywords['top'],
            bottom=keywords['bottom']
        )

        plt.savefig(filename + '.' + fileformat)
        plt.clf()
        log.file_written(filename + '.' + fileformat)

    def plot_two_concepts(
        self,
        concept,
        cogA,
        cogB,
        labels={1: '1', 2: '2', 3: '3', 4: '4'},
        tcolor={
            1: 'white',
            2: 'black',
            3: '0.5',
            4: '0.1'
        },
        filename='pdf',
        fileformat='pdf',
        threshold=1,
        usetex=True
    ):
        """
        Plot the evolution of two concepts in space.

        Notes
        -----
        This function may be useful to contrast patterns of different words in
        geographic space.

        """
        # usetex
        mpl.rc('text', usetex=True)

        # redefine taxa and tree for convenience
        taxa = self.taxa

        # XXX check for coordinates of the taxa, otherwise load them from file and
        # add them to the wordlist XXX add later, we first load it from file
        if 'coords' in self._meta:
            coords = self.coords

        else:
            coords = csv2dict(
                self.dataset,
                'coords',
                dtype=[str, float, float]
            )

        # load the rc-file XXX add internal loading later
        conf = self._config()

        # get the paps
        these_taxa = {}
        for taxon in taxa:

            # get the dictionary and the entry
            try:
                cogs = self.get_dict(col=taxon, entry='pap')[concept]
            except:
                cogs = []

            # check for identical cogs and assign them to the 4 categories
            if cogA in cogs and cogB in cogs:
                these_taxa[taxon] = 3
            elif cogA in cogs and cogB not in cogs:
                these_taxa[taxon] = 1
            elif cogA not in cogs and cogB in cogs:
                these_taxa[taxon] = 2
            else:
                these_taxa[taxon] = 4

        # determine the maxima of the coordinates
        latitudes = [i[0] for i in coords.values()]
        longitudes = [i[1] for i in coords.values()]

        min_lat, max_lat = min(latitudes), max(latitudes)
        min_lon, max_lon = min(longitudes), max(longitudes)

        # start to initialize the basemap
        fig = plt.figure()
        fig.add_subplot(111)

        # instantiate the basemap
        m = bmp.Basemap(
            llcrnrlon=min_lon + conf['min_lon'],
            llcrnrlat=min_lat + conf['min_lat'],
            urcrnrlon=max_lon + conf['max_lon'],
            urcrnrlat=max_lat + conf['max_lat'],
            resolution=conf['resolution'],
            projection=conf['projection']
        )

        # draw first values
        m.drawmapboundary(fill_color=conf['water_color'])
        m.drawcoastlines(color=conf['continent_color'], linewidth=0.5)
        m.drawcountries(color=conf['coastline_color'], linewidth=0.5)
        m.fillcontinents(color=conf['continent_color'], lake_color=conf['water_color'])

        # plot the points for the languages
        cell_text = []
        legend_check = []
        for i, (taxon, (lng, lat)) in enumerate(
                sorted(coords.items(), key=lambda x: x[0])):

            # retrieve x and y from the map
            x, y = m(lat, lng)

            # get the color of the given taxon
            # taxon_color = colors[groups[taxon]]

            if these_taxa[taxon] == 4:
                marker = '*'
            else:
                marker = 's'

            # check for legend
            if labels[these_taxa[taxon]] in legend_check:
                # check for forth marker
                # plot the marker
                plt.plot(
                    x,
                    y,
                    marker,
                    markersize=conf['markersize'],
                    color=tcolor[these_taxa[taxon]],
                    # zorder = 50,
                )
            else:
                # plot the marker
                plt.plot(
                    x,
                    y,
                    marker,
                    markersize=conf['markersize'],
                    color=tcolor[these_taxa[taxon]],
                    # zorder = 52,
                    label=labels[these_taxa[taxon]]
                )
                legend_check.append(labels[these_taxa[taxon]])

            # add number to celltext
            if usetex:
                cell_text.append([str(i + 1), taxon.replace('_', r'\_')])
            else:
                cell_text.append([str(i + 1), taxon])

            # plot the text
            if tcolor[these_taxa[taxon]] == 'black':
                textcolor = 'white'
            else:
                textcolor = 'black'

            plt.text(
                x,
                y,
                str(i + 1),
                size=str(int(conf['markersize'] / 2)),
                label=taxon,
                color=textcolor,
                horizontalalignment='center',
                verticalalignment='center',
            )

        this_table = plt.table(
            cellText=cell_text,
            colWidths=conf['table.column.width'],
            loc=conf['table.location'],
        )

        # adjust the table
        for line in this_table._cells:
            this_table._cells[line]._text._horizontalalignment = 'left'
            this_table._cells[line]._text._fontproperties.set_weight('bold')
            this_table._cells[line]._text.set_color(conf['table.text.color'])
            this_table._cells[line].set_height(conf['table.cell.height'])
            this_table._cells[line]._text._fontproperties.set_size(
                conf['table.text.size'])
            this_table._cells[line].set_linewidth(0.0)
            this_table._cells[line].set_color(conf['table.cell.color'])

        this_table.set_zorder(100)

        plt.legend(
            loc=conf['legend.location'],
            numpoints=1,
            prop={
                'size': conf['legend.size'],
                'weight': 'bold'
            }
        )

        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

        plt.savefig(filename + '.' + fileformat)
        plt.clf()

    def plot_GLS(
        self,
        glm,
        **keywords
    ):
        """
        Plot the inferred scenarios for a given model.
        """
        kw = dict(fileformat='png')
        kw.update(keywords)

        # store the graph
        for cog in self.cogs:
            gls = self.gls[glm][cog][0]
            g = gls2gml(
                gls,
                self.tgraph,
                self.tree,
                filename=self._output_path(
                    'gml', '{0}-{1}'.format(self.dataset, glm), cog))

            # if plot of gml is chosen
            nodes = []

            for n, d in g.nodes(data=True):
                x = d['graphics']['x']
                y = d['graphics']['y']
                f = d['graphics']['fill']
                try:
                    r = d['graphics']['angle']
                    s = d['graphics']['s']
                except:
                    r = None
                    s = None

                o = d['origin']
                l = d['label']

                nodes.append((x, y, f, o, l, r, s))

            edges = []
            for a, b, d in g.edges(data=True):
                xA = g.node[a]['graphics']['x']
                xB = g.node[b]['graphics']['x']
                yA = g.node[a]['graphics']['y']
                yB = g.node[b]['graphics']['y']

                edges += [(xA, xB, yA, yB)]

            # mpl.rc('text',usetex=keywords['usetex'])
            fig = plt.figure()
            fig.add_subplot(111)
            plt.axes(frameon=False)
            plt.xticks([])
            plt.yticks([])

            plt.axis('equal')

            for xA, xB, yA, yB in edges:
                plt.plot(
                    [xA, xB],
                    [yA, yB],
                    '-',
                    color='black',
                    linewidth=5
                )
                plt.plot(
                    [xA, xB],
                    [yA, yB],
                    '-',
                    color='0.2',
                    linewidth=4
                )
            for x, y, f, o, l, r, s in nodes:

                if f == '#000000':
                    f = '#a3a3a3'
                    c = '#a3a3a3'
                else:
                    c = '#000000'

                if o == 1:
                    size = 20
                else:
                    size = 10
                if l.startswith('edge') or l.startswith('root'):
                    plt.plot(x, y, 'o', markersize=size, color=f)
                else:
                    if not r:
                        plt.text(
                            x,
                            y,
                            l,
                            horizontalalignment='center',
                            verticalalignment='center',
                            size=8,
                            fontweight='bold',
                            color=c,
                            backgroundcolor=f
                        )
                    else:
                        plt.text(
                            x,
                            y,
                            l,
                            ha=s,
                            va='baseline',
                            size=8,
                            fontweight='bold',
                            color=c,
                            rotation=r,
                            rotation_mode='anchor',
                            bbox=dict(
                                facecolor='white',
                                boxstyle='square,pad=0.25',
                                ec="none",
                                alpha=0.25
                            ),
                        )

            plt.savefig(
                self._output_path(
                    'gml',
                    '{0}-{1}-figures'.format(self.dataset, glm),
                    '{0}-{1}.'.format(self.pap2con[cog], cog) + kw['fileformat']
                )
            )
            plt.clf()

    def get_stats(
        self,
        glm,
        subset='',
        filename=''
    ):
        """
        Calculate basic statistics for a given gain-loss model.
        """
        if not subset:
            gains = [b for a, b in self.gls[glm].values()]
        else:
            gains = []
            for cog in self.cogs:

                # get the respective subset-item first
                item = self[[c[0] for c in self.etd[cog] if c != 0][0], subset[0]]

                # check whether subset is as specified
                if item in subset[1] or item == subset[1]:
                    gains += [self.gls[glm][cog][1]]

        noo = sum(gains) / len(gains)

        ppc = sum([1 for g in gains if g > 1]) / len(gains)

        log.info('Number of Origins: {0:.2f}'.format(noo))
        log.info('Percentage of Patchy Cognates: {0:.2f}'.format(ppc))

        if not filename:
            return noo, ppc
        self._write_file(
            filename,
            'Number of origins: {0:.2f}\nPercentage of patchy cogs {1:.2f}\n'.format(noo,
                                                                                     ppc))

    def plot_concept_evolution(
        self,
        glm,
        concept='',
        fileformat='png',
        **keywords
    ):
        """
        Plot the evolution of specific concepts along the reference tree.
        """
        # make defaults
        defaults = dict(
            figsize=(15, 15),
            left=0.05,
            top=0.95,
            bottom=0.05,
            right=0.95,
            colormap=mpl.cm.jet,
            edgewidth=5,
            radius=2.5,
            outer_radius=0.5,
            inner_radius=0.25,
            cognates='',
            ref='',
            usetex=False,
            latex_preamble=False,
            textsize=8,
            subset=[]
        )

        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        # check for the correct item
        if not concept:
            concepts = self.concepts
        else:
            concepts = [i for i in self.concepts if i == concept]

        # XXX customize later XXX
        colormap = keywords['colormap']

        # start with the analysis
        for concept in concepts:
            log.info("Plotting concept '{0}'...".format(concept))

            # switch backend, depending on whether tex is used or not
            backend = mpl.get_backend()
            if keywords['usetex'] and backend != 'pgf':
                plt.switch_backend('pgf')
                mpl.rcParams['text.latex.unicode'] = True
            elif not keywords['usetex'] and backend != 'TkAgg':
                plt.switch_backend('TkAgg')

            # check for preamble settings
            if keywords['latex_preamble']:
                mpl.rcParams['pgf.preamble'] = keywords['latex_preamble']

            # make a graph
            graph = nx.Graph()

            # get all paps that are no singletons
            paps = sorted(set([p for p in self.get_list(
                row=concept,
                flat=True,
                entry='pap'
            ) if p not in self.singletons]))

            if len(paps) <= 0:
                log.warn(
                    "No entries for concept {0} could be found, skipping the plot.".format(
                        concept))
            else:

                # get the number of paps in order to get the right colors
                cfunc = np.array(np.linspace(10, 256, len(paps)), dtype='int')
                colors = dict([(paps[i], mpl.colors.rgb2hex(colormap(cfunc[i]))) for i in
                               range(len(paps))])

                # get the wedges for the paps
                wedges = {}
                linsp = np.linspace(0, 360, len(paps) + 1)
                for i, pap in enumerate(paps):
                    theta1, theta2 = linsp[i], linsp[i + 1]
                    wedges[pap] = (theta1, theta2)

                legendEntriesA = []
                legendTextA = []

                # add stuff for the legend
                for pap in paps:
                    w = mpl.patches.Wedge(
                        (0, 0),
                        1,
                        wedges[pap][0],
                        wedges[pap][1],
                        facecolor=colors[pap],
                        zorder=1,
                        linewidth=2,
                        edgecolor='black'
                    )
                    legendEntriesA += [w]
                    if keywords['cognates']:
                        idx = [x[0] for x in self.etd[pap] if x != 0][0]
                        legendTextA += [self[idx, keywords['cognates']]]
                    else:
                        legendTextA += [pap]

                # second legend explains evolution
                legendEntriesB = []
                legendTextB = []
                p = mpl.patches.Wedge(
                    (0, 0),
                    1,
                    0,
                    360,
                    facecolor='0.5',
                    linewidth=2,
                    edgecolor='black',
                )
                legendEntriesB += [p]
                legendTextB += ['Loss Event']
                p, = plt.plot(0, 0, '--', color='black', linewidth=2)
                legendEntriesB += [p]
                legendTextB += ['Gain Event']

                # overwrite stuff
                plt.plot(0, 0, 'o', markersize=2, zorder=2, color='white')

                # iterate over the paps and append states to the graph
                for pap in paps:

                    # get the graph with the model
                    gls = self.gls[glm][pap][0]
                    g = gls2gml(
                        gls,
                        self.tgraph,
                        self.tree,
                        filename=''
                    )

                    # iterate over the graph
                    for n, d in g.nodes(data=True):

                        # add the node if necessary
                        if n not in graph:
                            graph.add_node(n)

                        # add a pap-dictionary if it's not already there
                        if 'pap' not in graph.node[n]:
                            graph.node[n]['pap'] = {}

                        # add data
                        graph.node[n]['pap'][pap] = d['state']

                # create the figure
                fig = plt.figure(figsize=keywords['figsize'])
                figsp = fig.add_subplot(111)
                figsp.axes.get_xaxis().set_visible(False)
                figsp.axes.get_yaxis().set_visible(False)

                plt.axis('equal')

                xvals = []
                yvals = []

                # iterate over edges first
                for nA, nB in g.edges():
                    gA = g.node[nA]['graphics']
                    gB = g.node[nB]['graphics']
                    xA, yA = gA['x'], gA['y']
                    xB, yB = gB['x'], gB['y']

                    plt.plot(
                        [xA, xB],
                        [yA, yB],
                        '-',
                        color='black',
                        linewidth=keywords['edgewidth']
                    )

                # now iterate over the nodes
                for n, d in graph.nodes(data=True):
                    cpaps = d['pap']
                    x, y = g.node[n]['graphics']['x'], g.node[n]['graphics']['y']

                    # get z-value which serves as zorder attribute
                    try:
                        z = 6 * len(self.tree.getConnectingEdges('root', n))
                    except:
                        z = 0

                    xvals += [x]
                    yvals += [y]

                    # plot the default marker
                    plt.plot(
                        x,
                        y,
                        'o',
                        markersize=5,
                        color='black',
                        zorder=50
                    )
                    # check for origins in cpaps
                    if 'O' in cpaps.values():
                        w = mpl.patches.Wedge(
                            (x, y),
                            keywords['radius'] + keywords['outer_radius'],
                            0,
                            360,
                            facecolor='white',
                            zorder=57 + z,
                            linewidth=2.5,
                            linestyle='dashed',
                        )
                        figsp.add_artist(w)
                    elif 'o' in cpaps.values():
                        w = mpl.patches.Wedge(
                            (x, y),
                            keywords['radius'] + keywords['outer_radius'],
                            0,
                            360,
                            facecolor='white',
                            zorder=56 + z,
                            linewidth=2.5,
                            linestyle='solid',
                        )
                        figsp.add_artist(w)

                    if 'L' in cpaps.values() and 'O' in cpaps.values():
                        w = mpl.patches.Wedge(
                            (x, y),
                            keywords['radius'] + keywords['outer_radius'],
                            0,
                            360,
                            facecolor='0.5',
                            zorder=58 + z,
                            linewidth=2.5,
                            edgecolor='black',
                            linestyle='dashed'
                        )
                        figsp.add_artist(w)

                    elif "L" in cpaps.values():
                        w = mpl.patches.Wedge(
                            (x, y),
                            keywords['radius'] + keywords['outer_radius'],
                            0,
                            360,
                            facecolor='0.5',
                            zorder=59 + z,
                            linewidth=2.5,
                            edgecolor='black',
                        )
                        figsp.add_artist(w)

                    # plot all wedges
                    for pap in cpaps:

                        theta1, theta2 = wedges[pap]
                        color = colors[pap]

                        # check for characteristics of this pap
                        if cpaps[pap] == 'L':

                            w = mpl.patches.Wedge(
                                (x, y),
                                keywords['radius'],
                                theta1,
                                theta2,
                                facecolor=color,
                                zorder=61 + z,
                                alpha=0.25,
                                linewidth=2,
                                edgecolor='black',
                                linestyle='dotted'
                            )
                            figsp.add_artist(w)

                        elif cpaps[pap] == 'o':

                            w = mpl.patches.Wedge(
                                (x, y),
                                keywords['radius'],
                                theta1,
                                theta2,
                                facecolor=color,
                                zorder=61 + z,
                                linewidth=2,
                                edgecolor='black'
                            )
                            figsp.add_artist(w)

                        elif cpaps[pap] == 'O':

                            w = mpl.patches.Wedge(
                                (x, y),
                                keywords['radius'],
                                theta1,
                                theta2,
                                facecolor=color,
                                zorder=61 + z,
                                linewidth=2,
                                edgecolor='black',
                                linestyle='dashed'
                            )
                            figsp.add_artist(w)

                    # add number for node
                    if n in self.taxa:
                        plt.text(
                            x,
                            y,
                            n,
                            size=keywords['textsize'],
                            verticalalignment='baseline',
                            backgroundcolor='white',
                            horizontalalignment='center',
                            fontweight='bold',
                            color='black',
                            bbox=dict(
                                facecolor='white',
                                boxstyle='square,pad=0.25',
                                ec="none",
                                alpha=1
                            ),
                            zorder=300
                        )

                plt.xlim((min(xvals or [0]) - 10, max(xvals or [0]) + 10))
                plt.ylim((min(yvals or [0]) - 10, max(yvals or [0]) + 10))

                legend1 = plt.legend(legendEntriesA, legendTextA, loc='upper right',
                                     numpoints=1)
                plt.legend(legendEntriesB, legendTextB, loc='lower right')
                figsp.add_artist(legend1)

                plt.subplots_adjust(
                    left=keywords['left'],
                    right=keywords['right'],
                    top=keywords['top'],
                    bottom=keywords['bottom']
                )

                plt.savefig(
                    self._output_path(
                        'items',
                        '{0}-{1}'.format(self.dataset, glm),
                        concept.replace('/', '_') + '.' + fileformat
                    )
                )
                plt.close()

        # return the graph
        return

    # add an alias for backwards compatibility


TreBor = PhyBo
