from lingpy.algorithm import cluster_util
from unittest import TestCase

class Tests(TestCase):
    
    
    def test_valid_cluster(self):
        clrA = [0, 1, 2, 2]
        clrB = [1, 1, 0, 2]
        assert cluster_util.valid_cluster(clrA)
        assert not cluster_util.valid_cluster(clrB)
    
    
    def test_generate_all_clusters(self):
        assert len(set(cluster_util.generate_all_clusters(2))) == 2
        assert len(set(cluster_util.generate_all_clusters(3))) == 4
    
    
    def test_order_cluster(self):

        assert cluster_util.order_cluster([1, 0]) == [0, 1]

    def test_mutate_cluster(self):

        assert len(cluster_util.mutate_cluster([0, 0, 1])) == 3

    def test_generate_random_clusters(self):

        a = cluster_util.generate_random_cluster(3, bias=False)
        b = cluster_util.generate_random_cluster(2, bias='lumper')
        c = cluster_util.generate_random_cluster(4, bias='splitter')
        
        assert a != b != c

