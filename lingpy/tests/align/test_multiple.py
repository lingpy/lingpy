# *-* coding: utf-8 *-*
"""
Testing multiple module.
"""
from __future__ import print_function, division, unicode_literals
from unittest import TestCase

from lingpy.align import Multiple, mult_align


class Tests(TestCase):
    def setUp(self):
        self.seqs = ['waldemar', 'woldemort', 'vladimir']
        self.msa = Multiple(self.seqs)

    def test_mult_align(self):
        self.assertEqual(
            ''.join(mult_align(['waldemar', 'woldemort', 'vladimir'])[0]), 'wal-demar-')

    def test___get__(self):
        self.msa.lib_align()
        assert self.msa[0]
        assert self.msa[0, 0]

    def test_prog_align(self):
        self.msa.prog_align()
        assert self.msa.alm_matrix[0] == list('wal-demar-')

    def test_lib_align(self):
        self.msa.lib_align()
        assert self.msa.alm_matrix[0] == list('w-aldemar-')

    def test_get_pid(self):
        self.msa.prog_align()
        pid = int(self.msa.get_pid() * 100)
        assert pid == 63

    def test_swap_check(self):
        self.msa.prog_align()
        assert self.msa.swap_check()

    def test_iterate_all_sequences(self):
        self.msa.prog_align()
        first = ''.join(self.msa.alm_matrix[0])

        self.msa.iterate_all_sequences()
        secnd = ''.join(self.msa.alm_matrix[0])

        assert first == secnd

    def test_iterate_orphans(self):
        self.msa.prog_align()
        first = ''.join(self.msa.alm_matrix[0])

        self.msa.iterate_orphans(0.5)
        secnd = ''.join(self.msa.alm_matrix[0])

        assert first == secnd

    def test_iterate_similar_gap_sites(self):
        self.msa.prog_align()
        first = ''.join(self.msa.alm_matrix[0])

        self.msa.iterate_similar_gap_sites()
        secnd = ''.join(self.msa.alm_matrix[0])

        assert first == secnd

    def test_iterate_clusters(self):
        self.msa.prog_align()
        first = ''.join(self.msa.alm_matrix[0])

        self.msa.iterate_clusters(0.5)
        secnd = ''.join(self.msa.alm_matrix[0])

        assert first == secnd

    def test_sum_of_pairs(self):
        self.msa.prog_align()
        assert 8 > self.msa.sum_of_pairs() > 7

    def test_get_pairwise_alignments(self):
        self.msa.prog_align()
        self.msa.get_pairwise_alignments()
        assert hasattr(self.msa, 'alignments')

    def test_get_peaks(self):
        self.msa.prog_align()
        assert self.msa.get_peaks()[2] == 10

    def test_get_local_peaks(self):
        self.msa.prog_align()
        self.msa.get_local_peaks()

        assert self.msa.local[0] == 0
