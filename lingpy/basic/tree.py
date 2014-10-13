# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-08-15 13:11
# modified : 2014-08-15 13:11
"""
Basic module for the handling of language trees.
"""

__author__="Johann-Mattis List, Taraka Rama"
__date__="2014-08-15"

from ..thirdparty.cogent import LoadTree,PhyloNode
from ..algorithm import TreeDist
import random

class Tree(PhyloNode):
    """
    Basic class for the handling of phylogenetic trees.

    Parameters
    ----------
    tree : {str file}
        A string or a file containing trees in NEWICK format.
        
    """

    def __init__(self, tree):
        
        # this is an absolutely nasty hack, but it helps us to maintain
        # lingpy-specific aspects of cogent's trees and allows us to include
        # them in our documentation
        if type(tree) == str:
            tmp = LoadTree(treestring=tree)
        else:
            tmp = LoadTree(tree)
        for key,val in tmp.__dict__.items():
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
            return TreeDist.grf(str(self), str(other), distance = 'rf')
        elif distance == 'symmetric':
            return self.compareByPartitions(other, debug=debug)
        else:
            raise ValueError("[!] The distance you defined is not available.")
        
    def getDistanceToRoot(self, node):

        subtree = self.getNodeMatchingName(node)
        parent = ''
        counter = 0
        while parent != self.Name:
            parent = subtree.Parent.Name
            subtree = self.getNodeMatchingName(parent)
            counter += 1

        return counter
    
    def getTransitionMatrix(self):
        """
        Return a matrix with two columns for each node in the tree.

        """
        matrix = []
        for name,node in self.getNodesDict().items():
            if node.Parent:
                matrix += [[name, node.Parent.Name]]

        return matrix
    
    def createTransitionMatrix(self, states, tmat, random=True):
        
        # get number of states
        unique_states = sorted(set(states))
        
        max_state = max(states)

        matrix = []
        counter = 0
        for i in range(self._edge_len):
            if tmat[i][0] in self.taxa:
                matrix += [states[counter]]
                counter += 1
            else:
                matrix += [[randint(0,max_state),randint(0,max_state)]]
        return matrix


