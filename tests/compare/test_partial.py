import pytest

import lingpy.algorithm.extra
from lingpy.compare.partial import Partial, _get_slices


@pytest.fixture
def part(test_data):
    return Partial(
        str(test_data / 'partial_cognates.tsv'), segments='segments', split_on_tones=True)


@pytest.fixture
def part2(test_data):
    return Partial(str(test_data / 'partial_cognates-scored.tsv'), segments='segments')


def test__get_slices():
    a = _get_slices(list('ba²te²'), split_on_tones=True)
    b = _get_slices(list('ba²te²'), split_on_tones=False)
    assert a[0][1] == 3
    assert b[0][1] == 6


def test_get_partial_scorer(part2):
    part2.get_partial_scorer(runs=10)


def test_get_partial_matrices(part):
    for method in ['upgma', 'single', 'complete', 'ward', 'mcl']:
        matrix = list(part._get_partial_matrices(cluster_method=method, concept="bird"))[0]
        assert isinstance(matrix[0][0], (float, int))

    if lingpy.algorithm.extra.igraph:  # pragma: no cover
        for concept, tracer, matrix in part._get_partial_matrices(cluster_method='infomap'):
            assert isinstance(concept, str)
            assert [x[0] for x in tracer]


def test_partial_cluster(part, part2):
    with pytest.raises(ValueError):
        part.partial_cluster(cluster_method='upgmu')
    part.partial_cluster(
        method='sca',
        threshold=0.45,
        split_on_tones=True,
        post_processing=False,
        cluster_method='infomap' if lingpy.algorithm.extra.igraph else 'upgma',
        ref='parts1'
    )

    part.partial_cluster(method='sca', threshold=0.45,
            post_processing=False, cluster_method='mcl', ref='parts2',
            split_on_tones=True)
    part.partial_cluster(method='sca', threshold=0.45,
            split_on_tones=True,
            post_processing=False,
            cluster_method='upgma', ref='parts3')

    part2.partial_cluster(method='lexstat', threshold=0.6,
            cluster_method='single', post_processing=True, imap_mode=False,
            split_on_tones=True, ref='parts4')

    # high threshold to trigger post-processing movement
    part.partial_cluster(method='sca', threshold=0.9,
            split_on_tones=True,
            cluster_method='single', post_processing=True, imap_mode=False,
            ref='parts5')

    assert part[9, 'parts3'][0] == part[10, 'parts3'][0]
    assert part2[8, 'parts4'][1] == part2[10, 'parts4'][1]


def test_add_cognate_ids(part):
    part.partial_cluster(method='sca', threshold=0.45,
            split_on_tones=True,
            cluster_method='upgma', ref='parts3')
    part.add_cognate_ids('parts3', 'cogs1', idtype='strict')
    part.add_cognate_ids('parts3', 'cogs2', idtype='loose')

    assert part[9, 'cogs1'] == part[10, 'cogs1']
    with pytest.raises(ValueError):
        part.add_cognate_ids('parts3', 'cogs1', idtype='dummy')
