from __future__ import unicode_literals
import os

from six import text_type
from mock import patch

from lingpy.tests.util import WithTempDir


class Tests(WithTempDir):
    @patch('lingpy.algorithm.clustering.nx', new=None)
    def test_fuzzy1(self):
        from lingpy.algorithm.clustering import fuzzy

        assert fuzzy(0.5, [[1, 0], [0, 1]], ['a', 'b']) is None

    @patch('lingpy.algorithm.clustering.nx')
    def test_fuzzy2(self, nx):
        from lingpy.algorithm.clustering import fuzzy

        for method in 'upgma simple complete'.split():
            for revert in [True, False]:
                fuzzy(0.5, [[1, 0], [0, 1]], ['a', 'b'], method=method, revert=revert)

    def test_matrix2tree(self):
        from lingpy.algorithm.clustering import matrix2tree

        newick = text_type(self.tmp_path('t'))
        matrix2tree([[1, 0], [0, 1]], ['a', 'b'], filename=newick)
        assert os.path.exists(newick + '.nwk')
        matrix2tree([[1, 0], [0, 1]], ['a', 'b'], tree_calc='upgma')

    def test_matrix2groups(self):
        from lingpy.algorithm.clustering import matrix2groups

        for method in 'upgma mcl simple complete'.split():
            matrix2groups(0.5, [[1, 0], [0, 1]], ['a', 'b'], cluster_method=method)

    def test_link_clustering(self):
        from lingpy.algorithm.clustering import link_clustering

        for mt in {"distances", "similarities", "weights"}:
            link_clustering(0.5, [[1, 0], [0, 1]], ['a', 'b'], matrix_type=mt)

    def test_partition_density(self):
        from lingpy.algorithm.clustering import partition_density

        partition_density([[1, 0], [0, 1]], 0.5)

    def test_best_threshold(self):
        from lingpy.algorithm.clustering import best_threshold

        best_threshold([[1, 0], [0, 1]])
