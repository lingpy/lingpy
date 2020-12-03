import pytest

from lingpy import LexStat
from lingpy.compare.partial import Partial
from lingpy.evaluate.acd import bcubes, partial_bcubes, pairs, diff, \
    random_cognates, extreme_cognates, npoint_ap


@pytest.fixture
def lex(test_data):
    return LexStat(str(test_data / 'KSL.qlc'))


@pytest.fixture
def part(test_data):
    p = Partial(str(test_data / 'partial_cognates.tsv'), segments='segments')
    p.add_entries('pid1', 'partial_cognate_sets', lambda x: x)
    p.add_entries('pid2', 'partialids2', lambda x: [int(y) for y in x.split(' ')])
    return p


def test_bcubes(lex):
    res = bcubes(lex, test='cogid', pprint=False)
    assert res == pytest.approx((1.0, 1.0, 1.0))
    _ = bcubes(lex, 'cogid', 'cogid', pprint=True, per_concept=True)


def test_partial_bcubes(part):
    res = partial_bcubes(part, 'pid1', 'pid2', pprint=False)
    assert [round(x, 2) for x in res] == [0.92, 0.98, 0.95]

    _ = partial_bcubes(part, 'pid1', 'pid2', pprint=True)


def test_pairs(lex):
    res = pairs(lex, test='cogid', pprint=False)
    assert res == pytest.approx((1.0, 1.0, 1.0))


def test_diff(lex, tmppath):
    res = diff(lex, test='cogid', tofile=False, pprint=False)
    assert res == (pytest.approx((1.0, 1.0, 1.0)), pytest.approx((1.0, 1.0, 1.0)))
    lex.add_entries('cugid', 'cogid', lambda x: x + 1 if x % 2 else x * x)

    fname = str(tmppath / 'test_acd')
    _ = diff(lex, gold='cogid', test='cogid', filename=fname, pprint=False)
    d2 = diff(lex, gold='cugid', test='cogid', filename=fname, pprint=False, tofile=False)
    _ = diff(lex, gold='cugid', test='cogid', filename=fname, pprint=False, tofile=True)
    assert d2[0] != 1


def test_random_cognates(lex):
    random_cognates(lex, ref='randomid')
    assert 'randomid' in lex.header


def test_extreme_cognates(lex):
    extreme_cognates(lex, ref="lumperid", bias='lumper')
    assert lex[1, 'lumperid'] == lex[2, 'lumperid']
    extreme_cognates(lex, ref='splitterid', bias='splitter')
    assert lex[1, 'splitterid'] != lex[2, 'splitterid']
    with pytest.raises(ValueError):
        extreme_cognates(lex, bias='')


def test_npoint_ap():
    scores = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    scores2 = scores[::-1]
    cognates1 = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    cognates2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    cognates3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    cognates4 = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    np1 = npoint_ap(scores, cognates1)
    np2 = npoint_ap(scores, cognates2)
    np3 = npoint_ap(scores, cognates3)
    np4 = npoint_ap(scores2, cognates1)
    np5 = npoint_ap(scores2, cognates1, reverse=True)
    np6 = npoint_ap(scores, cognates4)
    np7 = npoint_ap(scores2, cognates4)

    assert np2 == 1
    assert np3 == 0
    assert np4 == 0.5
    assert np5 == np1
    assert np6 == 1.0
    assert round(np7, 2) == 0.35
