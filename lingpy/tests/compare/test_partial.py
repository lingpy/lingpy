from __future__ import print_function, division, unicode_literals
from lingpy.compare.partial import Partial
from nose.tools import assert_raises
from lingpy.tests.util import test_data, WithTempDir
from six import text_type
import lingpy.algorithm.extra


class Tests(WithTempDir):

    def setUp(self):
        WithTempDir.setUp(self)
        self.part = Partial(test_data('partial_cognates.tsv'), segments='segments')
        self.part2 = Partial(test_data('partial_cognates-scored.tsv'),
                segments='segments')

    def test_get_partial_matrices(self):

        for method in ['upgma', 'single', 'complete', 'ward', 'mcl']:
            matrix = list(self.part._get_partial_matrices(cluster_method=method,
                concept="bird"))[0]
            assert isinstance(matrix[0][0], (float, int))
        
        if lingpy.algorithm.extra.igraph:
            for concept, tracer, matrix in self.part._get_partial_matrices(
                    cluster_method='infomap'):
                assert isinstance(concept, text_type)
                assert [x[0] for x in tracer]


    def test_partial_cluster(self):
        
        assert_raises(ValueError, self.part.partial_cluster, cluster_method='upgmu')
        self.part.partial_cluster(
                method='sca', threshold=0.45, cluster_method='infomap'\
                        if lingpy.algorithm.extra.igraph else 'upgma',
                        ref='parts1')
        self.part.partial_cluster(
                method='sca', threshold=0.45, cluster_method='mcl',
                        ref='parts2')
        self.part.partial_cluster(
                method='sca', threshold=0.45, cluster_method='upgma',
                        ref='parts3')

        self.part2.partial_cluster(
                method='lexstat', threshold=0.6, cluster_method='single',
                post_processing=True, imap_mode=False, ref='parts4')
        # high threshold to trigger post-processing movement
        self.part.partial_cluster(
                method='sca', threshold=0.9, cluster_method='single',
                post_processing=True, imap_mode=False, ref='parts5')

        assert self.part[9, 'parts3'][0] == self.part[10, 'parts3'][0]
        assert self.part2[8, 'parts4'][1] == self.part2[10, 'parts4'][1]

    def test_add_cognate_ids(self):
        self.part.partial_cluster(
                method='sca', threshold=0.45, cluster_method='upgma',
                        ref='parts3')
        self.part.add_cognate_ids('parts3', 'cogs1', idtype='strict')
        self.part.add_cognate_ids('parts3', 'cogs2', idtype='loose')
        assert self.part[9,'cogs1'] == self.part[10, 'cogs1']
        assert_raises(ValueError, self.part.add_cognate_ids, 'parts3', 'cogs1', idtype='dummy')

