# *-* coding: utf-8 *-*
"""
Utility functions for borrowing detection with the PhyBo class.
"""
from __future__ import unicode_literals, division, print_function

from lingpy import log
from lingpy.util import write_text_file, as_string

try:
    import scipy.stats as sps
except ImportError:
    log.missing_module('scipy')


def tstats(
    wordlist,
    glm='',
    network=False,
    acs=False,
    tree=False,
    singletons=True,
    return_dists=False
):
    """
    Calculate transmission statistics for a given MLN.
    """

    # check for attributes
    # return if no glm and no network is given
    if not glm and not network:
        raise ValueError("You must specify at least one network or a gain-loss model.")

    # check for acs and network
    if glm:
        # network = wordlist.graph[glm]
        acs = wordlist.acs[glm]

    # check for tree
    if not tree:
        tree = wordlist.tree

    # add the distributions of the leaves to the acs
    for t in tree.taxa:
        paps = wordlist.get_list(taxa=t, entry='pap', flat=True)
        cons = wordlist.get_list(taxa=t, entry='concept', flat=True)

        acs[t] = [(p, c) for p, c in zip(paps, cons)]

    # now we apply a simple way to resolve directions by taking the first
    # occurence of links in the tree to be the innovation, and all dependent
    # links to be the source of borrowings

    # create a queue
    queue = ['root']

    # make dictionary of innovated chars: these are currently all present in the
    # root, we order list as [inheritance,innovation,transfer]
    tracer = dict([(c[0], [0, 1, 0]) for c in acs['root']])
    states = {}

    # start to iterate
    while queue:

        # get current node
        node = queue.pop(0)

        # get the children
        children = tree.getNodeMatchingName(node).Children

        # get the chars of the node
        node_chars = list(set([c[0] for c in acs[node]]))

        # if there are children
        for child in children:

            # get the node name
            name = child.Name

            # append name to the queue
            queue += [name]

            # get the chars of the child
            chars = list(set([c[0] for c in acs[name]]))

            inn = 0
            ret = 0
            bor = 0
            # iterate over chars and decide where they come from
            for char in chars:

                if char not in wordlist.singletons or not singletons:
                    # if char is inherited, increase the score
                    if char in node_chars:
                        tracer[char][0] += 1
                        ret += 1

                    # if occurs the first time, it is an innovation
                    elif char not in tracer:
                        tracer[char] = [0, 1, 0]
                        inn += 1

                    # if it is in the tracer
                    elif char not in node_chars and char in tracer:
                        tracer[char][2] += 1
                        bor += 1

            states[name] = [ret, inn, bor]

    # calculate the scores
    ret = sum([c[0] for c in tracer.values()])
    inn = sum([c[1] for c in tracer.values()])
    tra = sum([c[2] for c in tracer.values()])

    ipn = inn / len(acs)
    tpn = tra / len(acs)

    total2 = ipn + tpn

    log.info("Innovations: {0}, {1:.2f}, {2:.2f}".format(inn, ipn, ipn / total2))
    log.info("Transferred: {0}, {1:.2f}, {2:.2f}".format(tra, tpn, tpn / total2))

    if return_dists:
        leaves = []
        nodes = []
        for node in [n for n in tree.getNodeNames() if n != 'root']:
            innovations = states[node][1] + states[node][2]
            if node in tree.taxa:
                leaves += [innovations]
            else:
                nodes += [innovations]

        # evaluate using mwu
        p, r = sps.mstats.kruskalwallis(leaves, nodes)

        return p, r

    return inn, tra, tracer


def get_acs(wordlist, glm, homoplasy=0, **keywords):
    """
    Calculate ancestral cognate distributions.
    """
    gls = wordlist.gls[glm]

    queue = ['root']
    acs = dict(
        root=[]
    )
    for char in gls:
        positives = [x[0] for x in gls[char][0] if x[1] == 1]
        if 'root' in positives:
            acs['root'] += [char]

    while queue:

        parent = queue.pop(0)

        children = wordlist.tree.getNodeMatchingName(parent).Children

        for child in children:
            child_name = child.Name
            queue += [child_name]
            acs[child_name] = []

            for char in gls:

                positives = [x[0] for x in gls[char][0] if x[1] == 1]
                negatives = [x[0] for x in gls[char][0] if x[1] == 0]

                if child_name in positives:
                    try:
                        acs[child_name] += [char]
                    except KeyError:
                        acs[child_name] = [char]
                elif child_name in negatives:
                    pass
                elif char in acs[parent]:
                    try:
                        acs[child_name] += [char]
                    except KeyError:
                        acs[child_name] = [char]

    dist = []
    for t in acs:
        if t not in wordlist.taxa:
            acs[t] = sorted(set(acs[t]))
            forms = len(acs[t])
            dist += [int(forms - forms * homoplasy + 0.5)]

    return acs, dist


def check_stats(models, wordlist, filename='results.txt', pprint=False):
    results = []
    for m in models:
        p, z = tstats(wordlist, m, return_dists=True)
        results += [[m, p, z]]


    txt = ''
    for a, b, c in results:
        txt += '{0}\t{1:.2f}\t{2:.2f}\n'.format(a, b, c)
    as_string(txt, pprint)
    if filename: write_text_file(filename, txt)
