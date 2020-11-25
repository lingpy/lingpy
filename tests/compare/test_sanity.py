import pytest

from lingpy.basic.wordlist import Wordlist
from lingpy.compare import sanity as sn


@pytest.fixture
def wl(test_data):
    return Wordlist(str(test_data / 'KSL5.qlc'))

def test__mutual_coverage(wl):
    assert len(sn._mutual_coverage('Albanian', 'English', wl, 'concept')) == 2


def test__get_concepts(wl):
    assert len(sn._get_concepts(wl, 'concept')) == 6


def test_mutual_coverage(wl):
    assert len(sn.mutual_coverage(wl)['French']['Albanian']) == 3


def test_mutual_coverage_check(wl):
    assert not sn.mutual_coverage_check(wl, 3)


def test_mutual_coverage_subset(wl):
    a, b = sn.mutual_coverage_subset(wl, 3, concepts='concept')
    assert a == 3
    assert b[0][0] == 3
    assert b[0][1][0] == 'Albanian'


def test_synonymy(wl):
    syns = sn.synonymy(wl)
    assert max(syns.values()) == 1
