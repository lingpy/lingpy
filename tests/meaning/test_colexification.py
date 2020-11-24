"""
Tests for colexification module.
"""
import pytest

from lingpy.meaning.colexification import *
import lingpy.meaning.colexification as colx
from lingpy.basic.wordlist import Wordlist


@pytest.fixture
def wordlist(test_data):
    return Wordlist(str(test_data / 'colexification.tsv'))


@pytest.fixture
def cols(wordlist):
    return colx._get_colexifications(wordlist)


def test_colexification_network(test_data, tmppath):
    graph = colexification_network(Wordlist(str(test_data / 'colexification.tsv')))
    assert "hand" in graph and "arm" in graph

    graph = colexification_network(Wordlist(str(test_data / 'colexification.tsv')), bipartite=True)
    assert 'arm' in graph['l4.4'] and 'hand' in graph['l4.4']

    _ = colexification_network(
        Wordlist(str(test_data / 'colexification.tsv')),
        output="gml",
        filename=str(tmppath / "test"))


def test__get_colexifications(cols):
    assert len(cols[0]) == 5


def test__get_colexifications_by_taxa(cols):
    colt = colx._get_colexifications_by_taxa(cols)
    assert ('arm', 'hand') in colt['l1']


def test__get_statistics(wordlist):
    stats = colx._get_statistics(wordlist)
    assert stats['arm']['wordOcc'] == 6


def test__make_matrix(cols):
    colt = colx._get_colexifications_by_taxa(cols)
    taxa = ['l1', 'l2', 'l3', 'l4', 'l5', 'l6']

    matrix = colx._make_matrix(taxa, colt)
    assert matrix[0][0] == 0


def test__make_graph(cols):
    graph = colx._make_graph(cols, bipartite=True)
    assert 'arm' in graph['l4.4'] and 'hand' in graph['l4.4']


def test_compare_colexifications(wordlist):
    matrix = colx.compare_colexifications(wordlist)
    assert matrix[0][0] == 0


def test_evaluate_colexifications(cols, tmppath):
    graph = colx._make_graph(cols)
    _, _ = colx.evaluate_colexifications(graph, weight='wordWeight')
    colx.evaluate_colexifications(graph, weight='wordWeight', outfile=str(tmppath / 'test'))
