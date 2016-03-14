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
    return ''.join([t for t in taxon if t not in '()[]{},;:."' + "'"])


def _xml2dict(
    infile,
    tax_name='pri-name'
):
    """
    Convert xml-based MultiTree format to Newick format.

    Parameters
    ----------
    infile : str
        Name of the input file.
    tax_name : str (default="pri-name")
        Name of the value that shall be used to create the tree. Use "code" in
        order to get the iso-code, "pri-name" will get the primary name.

    Returns
    -------
    newick,taxa : dict,list
        A dictionary in tree structure and a list of the taxon-names.
    """

    # parse the xml-file
    document = {}

    # get the document
    document['document'] = minidom.parse(infile)

    # get the hash
    document['hash'] = document['document'].getElementsByTagName('hash')[0]

    # get the tree
    document['tree'] = document['hash'].getElementsByTagName('tree')[0]

    # get the root
    document['root'] = document['tree'].getElementsByTagName('root')[0]

    # now start iteration
    nwk = {(0, 'root'): []}
    queue = [(document['root'], (0, "root"))]
    taxa = []
    while queue:

        root, (idx, tname) = queue.pop()

        max_idx = max([k[0] for k in nwk if type(k) == tuple])  # type(k) == int])

        if (idx, tname) not in nwk:
            nwk[(idx, tname)] = []

        # get the children
        children = [c for c in root.childNodes if c.nodeName == 'children']

        # get the childs
        childs = [c for c in children[0].childNodes if c.nodeName == 'child']

        # print("Idx {0} has {1} childs".format(idx,len(childs)))

        if childs:
            # iterate over childs
            for i, child in enumerate(childs):

                # get the name of the child
                name = [c for c in child.childNodes if c.nodeName == tax_name]
                name = name[0].childNodes[0].data

                if idx < max_idx + i + 1:
                    queue += [(child, (max_idx + i + 1, name))]
                    nwk[idx, tname] += [(max_idx + i + 1, name)]
                else:
                    name = [c for c in child.childNodes if c.nodeName == tax_name]
                    name = name[0].childNodes[0].data
                    nwk[idx, tname] = [_nwk_format(name)]
                    taxa.append(_nwk_format(name))

        else:
            name = [c for c in root.childNodes if c.nodeName == tax_name][0]
            name = name.childNodes[0].data

            nwk[idx, tname] = [_nwk_format(name)]
            taxa.append(_nwk_format(name))

    return nwk, taxa


def xml2nwk(
    infile,
    filename='',
    tax_name='pri-name',
):
    """
    Convert xml-based MultiTree format to Newick-format.

    Parameters
    ----------
    infile : str
        Name of the input file.
    filename : str (default="")
        If string is empty, the data will be returned as string, if a full
        string is passed, the data will be written to a file with that name.
    tax_name : str (default="pri-name")
        Name of the value that shall be used to create the tree. Use "code" in
        order to get the iso-code, "pri-name" will get the primary name.

    Returns
    -------
    newick : str
        A newick-string, if filename is not left empty.

    """
    nwk, taxa = _xml2dict(infile, tax_name)

    # first, create a specific newick-dictionary
    newick = {}

    # make a lambda function for conversion of internal nodes into names
    makeChild = lambda x: '{{x_{0}_{1}}}'.format(x[0], _nwk_format(x[1])) if type(
        x) == tuple else x

    for i, n in sorted(nwk, key=lambda x: x[0]):  # range(len(nwk)):

        name = _nwk_format(n)
        # create format-string for children
        children = [makeChild(c) for c in nwk[i, n]]

        # create dictionary to replace previous format string
        if len(children) > 1:
            newick['x_' + str(i) + '_' + name] = '(' + ','.join(children) + ')' + name
        else:
            newick['x_' + str(i) + '_' + name] = children[0]

    # add the taxa
    for taxon in taxa:
        newick[str(taxon)] = taxon

    # create the newick-string
    newick_string = "{x_0_root};"
    newick_check = newick_string

    # start conversion
    i = 0
    while True:

        newick_string = newick_string.format(**newick)
        if newick_check == newick_string:
            break
        else:
            newick_check = newick_string

    return newick_string

    if not filename:
        return newick_string
    util.write_text_file(filename + '.nwk', newick_string)
    return


def nwk2guidetree(
    newick
):
    """
    Build a tree matrix for a guide tree given in Newick format.
    Input is a binary tree with zero-based integer names at the leaves.
    """
    # assumption: a binary tree with integer names starting with 0 at the leaves
    tree = cg.LoadTree(treestring=newick)
    nodeIndex = {}
    nextIdx = len(tree.tips())
    # generate virtual cluster IDs for the tree nodes, store them in nodeIndex
    for node in tree.postorder():
        if not node.isTip():
            nodeIndex[node] = nextIdx
            nextIdx += 1
        else:
            nodeIndex[node] = int(node.Name)
    # construct tree matrix by another postorder traversal
    tree_matrix = []
    queue = deque(tree.postorder())
    while len(queue) > 0:
        curNode = queue.popleft()
        if not curNode.isTip():
            leftChild = curNode.Children[0]
            rightChild = curNode.Children[1]
            tree_matrix.append([nodeIndex[leftChild], nodeIndex[rightChild], 0.5, 0.5])
    return tree_matrix


def selectNodes(tree, selIndices):
    selNodes = []
    for leaf in tree.tips():
        if int(leaf.Name) in selIndices:
            selNodes.append(leaf)
    return selNodes


def treePath(node):
    path = [node]
    while not node.isRoot():
        node = node.Parent
        path.insert(0, node)
    return path


def constructSubtree(paths, index, curNode, indexMap):
    # create a map [node -> all paths containing that node at index position]
    partition = {}
    for node in {path[index] for path in paths}:
        partition[node] = [path for path in paths if path[index] == node]
    if len(partition) == 1:
        # no split, we simply go on to the next index in paths
        constructSubtree(paths, index + 1, curNode, indexMap)
    else:
        # split according to the partition, creating a new node where necessary
        for node in partition.keys():
            if len(partition[node]) == 1:
                # we have arrived at a leaf (or a unary branch above it), copy the leaf
                newLeafName = str(indexMap[int(partition[node][0][-1].Name)])
                newLeaf = cg.tree.TreeNode(Name=newLeafName)
                newLeaf.orig = partition[node][0][-1]
                curNode.Children.append(newLeaf)
                newLeaf.Parent = curNode
            else:
                newNode = cg.tree.TreeNode()
                newNode.orig = node
                curNode.Children.append(newNode)
                newNode.Parent = curNode
                constructSubtree(partition[node], index + 1, newNode, indexMap)


def subGuideTree(tree, selIndices):
    selNodes = selectNodes(tree, selIndices)
    indexMap = dict(zip(selIndices, range(0, len(selIndices))))
    paths = [treePath(node) for node in selNodes]
    # print str(paths)
    subtree = cg.tree.TreeNode()
    subtree.orig = tree.root()
    constructSubtree(paths, 1, subtree, indexMap)
    return subtree


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
