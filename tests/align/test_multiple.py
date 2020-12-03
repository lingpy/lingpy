"""
Testing multiple module.
"""
import pytest

from lingpy.align import Multiple, mult_align


@pytest.fixture
def seqs():
    return ['waldemar', 'woldemort', 'vladimir']


@pytest.fixture
def msa(seqs):
    return Multiple(seqs) 


def test_mult_align(seqs):
    assert ''.join(mult_align(seqs)[0]) == 'wal-demar-'


def test___get__(msa):
    msa.lib_align()
    assert msa[0]
    assert msa[0, 0]


def test_prog_align(msa):
    msa.prog_align()
    assert msa.alm_matrix[0] == list('wal-demar-')
    msa.prog_align(sonars=[
        [1, 2, 3, 4, 1, 2, 3, 4],
        [1, 2, 3, 4, 5, 1, 2, 3, 4],
        [1, 1, 1, 1, 1, 1, 1, 1]])
    assert msa.alm_matrix[0] == list('wal-demar-')


def test_lib_align(msa):
    msa.lib_align()
    assert msa.alm_matrix[0] == list('w-aldemar-')


def test_get_pid(msa):
    msa.prog_align()
    pid = int(msa.get_pid() * 100)
    assert pid == 63


def test_swap_check(msa):
    msa.prog_align()
    assert msa.swap_check()


def test_iterate_all_sequences(msa):
    msa.prog_align()
    first = ''.join(msa.alm_matrix[0])

    msa.iterate_all_sequences()
    secnd = ''.join(msa.alm_matrix[0])

    assert first == secnd


def test_iterate_orphans(msa):
    msa.prog_align()
    first = ''.join(msa.alm_matrix[0])

    msa.iterate_orphans(0.5)
    secnd = ''.join(msa.alm_matrix[0])

    assert first == secnd


def test_iterate_similar_gap_sites(msa):
    msa.prog_align()
    first = ''.join(msa.alm_matrix[0])

    msa.iterate_similar_gap_sites()
    secnd = ''.join(msa.alm_matrix[0])

    assert first == secnd


def test_iterate_clusters(msa):
    msa.prog_align()
    first = ''.join(msa.alm_matrix[0])

    msa.iterate_clusters(0.5)
    secnd = ''.join(msa.alm_matrix[0])

    assert first == secnd


def test_sum_of_pairs(msa):
    msa.prog_align()
    assert 8 > msa.sum_of_pairs() > 7


def test_get_pairwise_alignments(msa):
    msa.prog_align()
    msa.get_pairwise_alignments()
    assert hasattr(msa, 'alignments')


def test_get_peaks(msa):
    msa.prog_align()
    assert msa.get_peaks()[2] == 10


def test_get_local_peaks(msa):
    msa.prog_align()
    msa.get_local_peaks()

    assert msa.local[0] == 0
