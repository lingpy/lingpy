import pytest

from lingpy import LexStat, Wordlist
from lingpy.compare.partial import Partial
from lingpy.evaluate.acd import (
        bcubes, partial_bcubes, pairs, diff,
        random_cognates, extreme_cognates, npoint_ap, _get_bcubed_score
        )


@pytest.fixture
def lex(test_data):
    return LexStat(str(test_data / 'KSL.qlc'))


@pytest.fixture
def part(test_data):
    p = Partial(str(test_data / 'partial_cognates.tsv'), segments='segments')
    p.add_entries('pid1', 'partial_cognate_sets', lambda x: x)
    p.add_entries('pid2', 'partialids2', lambda x: [int(y) for y in x.split(' ')])
    return p


def test_get_bcubed_score():
    a, b = [1, 1, 2, 2], [1, 1, 1, 1]
    assert _get_bcubed_score(a, b) == 1


def test_bcubes(lex):
    res = bcubes(lex, test='cogid', pprint=False)
    assert res == pytest.approx((1.0, 1.0, 1.0))
    _ = bcubes(lex, 'cogid', 'cogid', pprint=True, per_concept=True)


@pytest.mark.parametrize('modify_ref', [
    lambda value: value + 10000,
    str,
    lambda value: value % 5 + 1,
])
@pytest.mark.parametrize('synonym', [False, True])
def test_bcubes_normalizes_cognates_without_changing_word_ids(modify_ref, synonym):
    wordlist = Wordlist({
        0: ['doculect', 'concept', 'ipa', 'cogid', 'testid'],
        1: ['A', 'hand', 'a', 1, 3],
        2: ['B', 'hand', 'b', 2, 1],
        3: ['C', 'hand', 'c', 2, 1],
        4: ['D', 'hand', 'd', 1, 1],
        5: ['E', 'hand', 'e', 1, 1],
        **({6: ['A', 'hand', 'a2', 1, 3]} if synonym else {}),
    })
    expected = (3 / 5, 11 / 15, 33 / 50)
    assert bcubes(wordlist, test='testid', pprint=False) == pytest.approx(expected)
    assert bcubes(wordlist, test='testid', modify_ref=modify_ref,
                  pprint=False) == pytest.approx(expected)


def test_partial_bcubes(part):
    res = partial_bcubes(part, 'pid1', 'pid2', pprint=False)
    assert [round(x, 2) for x in res] == [0.92, 0.98, 0.95]

    _ = partial_bcubes(part, 'pid1', 'pid2', pprint=True)


def test_pairs(lex):
    res = pairs(lex, test='cogid', pprint=False)
    assert res == pytest.approx((1.0, 1.0, 1.0))


def test_diff(lex, tmp_path):
    res = diff(lex, test='cogid', tofile=False, pprint=False)
    assert res == (pytest.approx((1.0, 1.0, 1.0)), pytest.approx((1.0, 1.0, 1.0)))
    lex.add_entries('cugid', 'cogid', lambda x: x + 1 if x % 2 else x * x)

    fname = str(tmp_path / 'test_acd')
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


@pytest.mark.parametrize('gold_ids', [(-1, 1, -1, 2), (1, -1, 1, -2)])
@pytest.mark.parametrize('per_concept', [False, True])
def test_bcubes_modify_loan_ids(gold_ids, per_concept):
    data = {0: ['doculect', 'concept', 'ipa', 'cogid', 'testid']}
    for i, (gold_id, test_id) in enumerate(zip(gold_ids, [4, 4, 5, 6]), 1):
        data[i] = ['language' + str(i), 'hand', 'a', gold_id, test_id]
    wordlist = Wordlist(data)
    assert bcubes(wordlist, test='testid', modify_ref=abs, pprint=False,
                  per_concept=per_concept) == pytest.approx((1., 2. / 3., .8))


@pytest.mark.parametrize('gold_ids', [(-1, 1, -1, 2), (1, -1, 1, -2)])
def test_diff_modify_loan_ids(gold_ids):
    data = {0: ['doculect', 'concept', 'ipa', 'cogid', 'testid']}
    for i, (gold_id, test_id) in enumerate(zip(gold_ids, [4, 4, 5, 6]), 1):
        data[i] = ['language' + str(i), 'hand', 'a', gold_id, test_id]
    wordlist = Wordlist(data)
    result = diff(wordlist, test='testid', modify_ref=abs, pprint=False, tofile=False)
    assert result[0] == pytest.approx((1., 2. / 3., .8))
    assert result[1] == pytest.approx((1., 1. / 3., .5))


@pytest.mark.parametrize('bias', [False, 'lumper', 'splitter'])
def test_random_cognates_separates_concepts(bias):
    wordlist = Wordlist({
        0: ['doculect', 'concept', 'ipa'],
        1: ['English', 'hand', 'h a n d'],
        2: ['English', 'foot', 'f u t'],
        3: ['English', 'head', 'h e d'],
    })
    random_cognates(wordlist, bias=bias)
    assert len({wordlist[idx, 'randomid'] for idx in wordlist}) == 3
    assert all(wordlist[idx, 'randomid'] > 0 for idx in wordlist)
