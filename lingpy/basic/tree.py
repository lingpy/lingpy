"""
Basic module for the handling of language trees.
"""
from __future__ import unicode_literals, division, print_function
import random

from lingpy.thirdparty.cogent import LoadTree, PhyloNode
from lingpy.algorithm import TreeDist


def _star_tree(taxa_list):
    """
    Initialize a star tree for the random tree function.
    """
    return "(" + ",".join(taxa_list) + ");"


def random_tree(taxa, branch_lengths=False):
    """
    Create a random tree from a list of taxa.

    Parameters
    ----------


    taxa : list
        The list containing the names of the taxa from which the tree will be
        created.
    branch_lengths : bool (default=False)
        When set to *True*, a random tree with random branch lengths will be
        created with the branch lengths being in order of the maximum number of
        the total number of internal branches.

    Returns
    -------
    tree_string : str
        A string representation of the random tree in Newick format.

    """
    # clone the list in order to avoid that lists used outside the function
    # suffer from modifications
    taxa_list = [t for t in taxa]

    random.shuffle(taxa_list)

    if not branch_lengths:
        while(len(taxa_list) > 1):
            ulti_elem = str(taxa_list.pop())
            penulti_elem = str(taxa_list.pop())
            taxa_list.insert(0, "(" + penulti_elem + "," + ulti_elem + ")")
            random.shuffle(taxa_list)

        taxa_list.append(";")
        return "".join(taxa_list)

    else:
        brlen_taxa_list = []
        nbr = 2 * len(taxa_list) - 3
        for taxon in taxa_list:
            brlen_taxa_list.append(str(taxon) + ":" + '{0:.2f}'.format(
                random.uniform(1, nbr)))
        while(len(brlen_taxa_list) > 1):
            ulti_elem = str(brlen_taxa_list.pop())
            penulti_elem = str(brlen_taxa_list.pop())
            if len(brlen_taxa_list) > 0:
                brlen_taxa_list.insert(
                    0,
                    "(" + penulti_elem + "," + ulti_elem + ")" + ":" + '{0:.2f}'.format(
                        random.uniform(0, nbr)))
            else:
                brlen_taxa_list.insert(0, "(" + penulti_elem + "," + ulti_elem + ")")
            random.shuffle(brlen_taxa_list)
        brlen_taxa_list.append(";")
        return "".join(brlen_taxa_list)


class Tree(PhyloNode):
    """
    Basic class for the handling of phylogenetic trees.

    Parameters
    ----------
    tree : {str file list}
        A string or a file containing trees in Newick format. As an
        alternative, you can also simply pass a list containing taxon names. In
        that case, a random tree will be created from the list of taxa.


    branch_lengths : bool (default=False)
        When set to *True*, and a list of taxa is passed instead of a Newick
        string or a file containing a Newick string, a random tree with random
        branch lengths will be created with the branch lengths being in order
        of the maximum number of the total number of internal branches.


    """

    def __init__(self, tree, **keywords):

        # this is an absolutely nasty hack, but it helps us to maintain
        # lingpy-specific aspects of cogent's trees and allows us to include
        # them in our documentation
        if type(tree) == str:
            if tree[-4:] not in ['.nwk', '.txt']:
                tmp = LoadTree(treestring=tree)
            else:
                tmp = LoadTree(tree)
        else:
            if type(tree) == list:
                tmp = LoadTree(treestring=random_tree(tree, **keywords))
            else:
                tmp = LoadTree(tree)

        for key, val in tmp.__dict__.items():
            self.__dict__[key] = val

        self._edge_len = len(self.getNodeNames()) - 1

    def get_distance(self, other, distance='grf', debug=False):
        """
        Function returns the Robinson Fould distance between the two trees.

        Parameters
        ----------
        other : lingpy.basic.tree.Tree
            A tree object. It should have the same number of taxa as the
            intitial tree.
        distance : { "grf", "rf", "branch", "symmetric"} (default="grf")
            The distance which shall be calculated. Select between:

            * "grf": the generalized Robinson Fould Distance
            * "rf": the Robinson Fould Distance
            * "branch": the distance in terms of branch lengths
            * "symmetric": the symmetric difference between all partitions of
                the trees


        """
        if distance == 'grf':
            return TreeDist.grf(str(self), str(other), distance='grf')
        elif distance in ['branch', 'branch_length', 'branchlength']:
            return self.distance(other)
        elif distance == 'rf':
            return TreeDist.grf(str(self), str(other), distance='rf')
        elif distance == 'symmetric':
            return self.compareByPartitions(other, debug=debug)
        else:
            raise ValueError("[!] The distance you defined is not available.")

    def getDistanceToRoot(self, node):
        """
        Return the distance from the given node to the root.
        
        Parameters
        ----------
        node : str
            The name of a given node in the tree.

        Returns
        -------
        distance : int
            The distance of the given node to the root of the tree.
        """

        subtree = self.getNodeMatchingName(node)
        parent = ''
        counter = 0
        while parent != self.Name:
            parent = subtree.Parent.Name
            subtree = self.getNodeMatchingName(parent)
            counter += 1

        return counter

    #def getTransitionMatrix(self):
    #    """
    #    Return a matrix with two columns for each node in the tree.

    #    """
    #    matrix = []
    #    for name, node in self.getNodesDict().items():
    #        if node.Parent:
    #            matrix += [[name, node.Parent.Name]]

    #    return matrix

    #def createTransitionMatrix(self, states, tmat):
    #    max_state = max(states)

    #    matrix = []
    #    counter = 0
    #    for i in range(self._edge_len):
    #        if tmat[i][0] in self.taxa:
    #            matrix += [states[counter]]
    #            counter += 1
    #        else:
    #            matrix += [[random.randint(0, max_state), random.randint(0, max_state)]]
    #    return matrix
