from __future__ import unicode_literals

import os

import pytest

from lingpy.algorithm.clustering import best_threshold, check_taxon_names, \
    find_threshold, flat_cluster, link_clustering, matrix2groups, matrix2tree, \
    neighbor, partition_density, upgma, fuzzy


@pytest.fixture(name='matrix')
def fixture_matrix():
    return [[0.0, 0.5, 0.67, 0.8, 0.2],
            [0.5, 0.0, 0.4, 0.7, 0.6],
            [0.67, 0.4, 0.0, 0.8, 0.8],
            [0.8, 0.7, 0.8, 0.0, 0.3],
            [0.2, 0.6, 0.8, 0.3, 0.0]]


@pytest.fixture(name='taxa')
def fixture_taxa():
    return ['German', 'Swedish', 'Icelandic', 'English', 'Dutch']


def test_check_taxa(taxa):
    with pytest.raises(ValueError):
        check_taxon_names(['Eng:lish', 'Ch!nese', 'G3rm4n'])

    try:
        check_taxon_names(taxa)
    except ValueError:
        pytest.fail()


def test_upgma(matrix, taxa):
    tree = upgma(matrix, taxa, distances=True)
    assert 'English' in tree

    with pytest.raises(ValueError):
        _ = upgma([[0, 1], [1, 0]], ['Eng:lish', 'Ger)man'])


def test_neighbor(matrix, taxa):
    tree = neighbor(matrix, taxa, distances=True)
    assert 'English' in tree

    with pytest.raises(ValueError):
        _ = neighbor([[0, 1], [1, 0]], ['Eng:lish', 'Ger)man'])


def test_fuzzy(matrix, taxa):
    for method in 'upgma simple complete'.split():
        for revert in [True, False]:
            _ = fuzzy(0.5, matrix, taxa, method=method, revert=revert)


def test_matrix2tree(matrix, taxa, tmpdir):
    newick = tmpdir.join('newick')
    _ = matrix2tree(matrix, taxa, filename=newick)
    assert os.path.exists(newick + '.nwk')

    _ = matrix2tree(matrix, taxa, tree_calc='upgma')
    _ = matrix2tree(matrix, taxa, tree_calc='neighbor')

    with pytest.raises(ValueError):
        _ = matrix2tree(matrix, taxa, tree_calc='does_not_exist')

    with pytest.raises(ValueError):
        _ = matrix2tree(*[matrix, taxa], **{"tree_calc": "dummy"})


def test_matrix2groups(matrix, taxa):
    for method in 'upgma mcl simple complete'.split():
        _ = matrix2groups(0.5, matrix, taxa, cluster_method=method)


def test_link_clustering(matrix, taxa):
    similarity_matrix = [[1 - cell for cell in row] for row in matrix]

    for r in [True, False]:
        for f in [True, False]:
            _ = link_clustering(0.5, matrix, taxa, matrix_type="distances",
                                revert=r, fuzzy=f)
            _ = link_clustering(0.5, similarity_matrix, taxa,
                                matrix_type="similarities", revert=r, fuzzy=f)
            _ = link_clustering(0.5, similarity_matrix, taxa,
                                matrix_type="weights", revert=r, fuzzy=f)

    with pytest.raises(ValueError):
        _ = link_clustering(0.5, matrix, taxa, matrix_type="dummy")


def test_partition_density(matrix):
    _ = partition_density(matrix, 0.5)


def test_best_threshold(matrix):
    _ = best_threshold(matrix, trange=(0.0, 1.0, 0.05))


def test_find_threshold(matrix):
    _ = find_threshold(matrix)
    _ = find_threshold(matrix, logs=False)
    assert find_threshold([[0, 1], [1, 0]]) is None


def test_check_taxon_names():
    with pytest.raises(ValueError):
        check_taxon_names(['k,k'])


def test_flat_cluster(matrix, taxa):
    for method in ['upgma', 'single', 'complete', 'ward']:
        _ = flat_cluster(method, 0.5, matrix, taxa, revert=True)
        _ = flat_cluster(method, 0.5, matrix, taxa, revert=False)
        _ = flat_cluster(method, 0.5, matrix, False, revert=False)
