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
        tmp = LoadTree(tree)
        for key,val in tmp.__dict__.items():
            self.__dict__[key] = val

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
    

