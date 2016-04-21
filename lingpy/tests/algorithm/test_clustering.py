from __future__ import unicode_literals
import os

from six import text_type
from mock import patch
from nose.tools import assert_raises

from lingpy.tests.util import WithTempDir

class Tests(WithTempDir):

    def setUp(self):
        WithTempDir.setUp(self)
        
        self.matrix = [[0.0, 0.5, 0.67, 0.8, 0.2],
                [0.5, 0.0, 0.4, 0.7, 0.6],
                [0.67, 0.4, 0.0, 0.8, 0.8],
                [0.8, 0.7, 0.8, 0.0, 0.3],
                [0.2, 0.6, 0.8, 0.3, 0.0]]
        self.taxa = ['German','Swedish','Icelandic','English','Dutch']

    def test_fuzzy(self):
        from lingpy.algorithm.clustering import fuzzy
        for method in 'upgma simple complete'.split():
            for revert in [True, False]:
                fuzzy(0.5, self.matrix, self.taxa, method=method, revert=revert)

    def test_matrix2tree(self):
        from lingpy.algorithm.clustering import matrix2tree

        newick = text_type(self.tmp_path('t'))
        matrix2tree(self.matrix, self.taxa, filename=newick)
        assert os.path.exists(newick + '.nwk')
        matrix2tree(self.matrix, self.taxa, tree_calc='upgma')
        matrix2tree(self.matrix, self.taxa, tree_calc='neighbor')

        assert_raises(ValueError, matrix2tree, *[self.matrix, self.taxa],
                **{"tree_calc": "dummy"})

    def test_matrix2groups(self):
        from lingpy.algorithm.clustering import matrix2groups

        for method in 'upgma mcl simple complete'.split():
            matrix2groups(0.5, self.matrix, self.taxa, cluster_method=method)

    def test_link_clustering(self):
        from lingpy.algorithm.clustering import link_clustering

        similarity_matrix = [[1-cell for cell in row] for row in self.matrix]
        
        for val in [True,False]:
            for val2 in [True, False]:
                link_clustering(0.5, self.matrix, self.taxa,
                        matrix_type="distances", revert=val, fuzzy=val2)
                link_clustering(0.5, similarity_matrix, self.taxa,
                        matrix_type="similarities", revert=val, fuzzy=val2)
                link_clustering(0.5, similarity_matrix, self.taxa,
                        matrix_type="weights", revert=val, fuzzy=val2)

        assert_raises(ValueError, link_clustering, 0.5, self.matrix, self.taxa,
                matrix_type="dummy")

    def test_partition_density(self):
        from lingpy.algorithm.clustering import partition_density

        partition_density(self.matrix, 0.5)

    def test_best_threshold(self):
        from lingpy.algorithm.clustering import best_threshold

        best_threshold(self.matrix, trange=(0.0, 1.0, 0.05))

    def test_find_threshold(self):

        from lingpy.algorithm.clustering import find_threshold
        find_threshold(self.matrix)
        find_threshold(self.matrix, logs=False)
        assert find_threshold([[0,1],[1,0]]) is None

    def test_flat_cluster(self):
        from lingpy.algorithm.clustering import flat_cluster
        for method in ['upgma', 'single', 'complete', 'ward']:
            flat_cluster(method, 0.5, self.matrix, self.taxa, revert=True)
            flat_cluster(method, 0.5, self.matrix, self.taxa, revert=False)
            flat_cluster(method, 0.5, self.matrix, False, revert=False)



