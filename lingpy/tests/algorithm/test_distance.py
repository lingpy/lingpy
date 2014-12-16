from __future__ import print_function, unicode_literals
import math

from unittest import TestCase


class Tests(TestCase):
    def test_hamming(self):
        from lingpy.algorithm.distance import hamming

        l1 = range(5)
        self.assertEqual(hamming(l1, l1), 0)
        self.assertEqual(hamming(l1, reversed(l1)), 4)
        self.assertEqual(hamming(l1, []), 5)

    def test_jaccard(self):
        from lingpy.algorithm.distance import jaccard

        self.assertEqual(jaccard(set(), set()), 0)
        jaccard({1}, {2})

    def test_euclidean(self):
        from lingpy.algorithm.distance import euclidean

        self.assertAlmostEqual(euclidean((0, 1), (1, 0)), math.sqrt(2))

    def test_pearson(self):
        from lingpy.algorithm.distance import pearson

        self.assertRaises(NotImplementedError, pearson, 1, 2)