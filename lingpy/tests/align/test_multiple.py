# *-* coding: utf-8 *-*
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-12-02 15:37
# modified : 2014-12-02 15:37
"""
Testing multiple module.
"""

__author__="Johann-Mattis List"
__date__="2014-12-02"

import os
import unittest
from lingpy.align import Multiple, mult_align

def test_mult_align():

    m = mult_align(['waldemar','woldemort','vladimir'])

    assert ''.join(m[0]) == 'wal-demar-'
    
class TestMultiple(object):

    def setup(self):

        self.seqs = ['waldemar','woldemort','vladimir']
        self.msa = Multiple(self.seqs)

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

        assert hasattr(self.msa,'alignments')

    def test_get_peaks(self):

        self.msa.prog_align()
        assert self.msa.get_peaks()[2] == 10

    def test_get_local_peaks(self):

        self.msa.prog_align()
        self.msa.get_local_peaks()

        assert self.msa.local[0] == 0


