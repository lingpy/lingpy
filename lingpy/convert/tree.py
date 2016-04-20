"""
Functions for tree calculations and working with trees.
"""
import xml.dom.minidom as minidom
from collections import deque

from lingpy.thirdparty import cogent as cg
from lingpy import util


def _nwk_format(taxon):
    """
    Function cleans a taxon string in order to make it apt for newick-representation.
    """
    # strip whitespace off
    taxon = taxon.strip()

    # replace white space underscore
    taxon = taxon.replace(' ', '_')

    # exclude all kinds of brackets
    return ''.join([t for t in taxon if t not in '!?()[]{},;:."' + "'"])

def nwk2tree_matrix(newick):
    """
    Convert a newick file to a tree matrix.

    Notes
    -----
    This is an additional function that can be used for plots with help of
    matplotlibs functions. The tree_matrix is compatible with those matrices
    that scipy's linkage functions create.
    """
    if type(newick) == str:
        tree = cg.LoadTree(treestring=newick)
    elif hasattr(newick, 'root'):
        tree = newick

    taxa = [t for t in sorted(
        tree.taxa,
        key=lambda x: len(tree.getConnectingEdges('root', x)),
        reverse=True
    )]

    tax2id = dict(zip(taxa, range(len(taxa))))
    nodes = [t for t in tree.getNodeNames() if t not in taxa]

    nodes = sorted(
        nodes,
        key=lambda x: len(tree.getNodeMatchingName(x).tips()),
    )
    matrix = []

    for node in nodes:
        n = tree.getNodeMatchingName(node)
        children = n.Children
        names = [c.Name for c in children]
        idxA = tax2id[names[0]]
        idxB = tax2id[names[1]]
        idx = max(tax2id.values()) + 1
        tax2id[node] = idx
        obs = len(n.tips())
        dst = obs * 1.0
        matrix += [[idxA, idxB, dst, obs]]

    return matrix, taxa
