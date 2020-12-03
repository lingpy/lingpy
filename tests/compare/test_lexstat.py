import os
import pathlib

import pytest
from clldutils import jsonlib

from lingpy import LexStat, rc
from lingpy.compare.lexstat import char_from_charstring, get_score_dict


def test_char_from_charstring():
    assert char_from_charstring('a.b.c') == "b"
    assert char_from_charstring('a.b') == "a"
    with pytest.raises(ValueError):
        char_from_charstring("a")


def test_get_score_dict():
    chars = ["1.A.-", "2.B.-"]
    model = rc("sca")
    sd = get_score_dict(chars, model)
    assert sd['A', 'B'] == -22.5


@pytest.fixture
def lextstat_factory(tmppath):
    def make(*args, **kw):
        kw.setdefault('errors', str(tmppath / 'errors.log'))
        return LexStat(*args, **kw)
    return make


@pytest.fixture
def lex(test_data, lextstat_factory):
    return lextstat_factory(str(test_data / 'KSL.qlc'))


@pytest.fixture
def log(mocker):
    return mocker.Mock(
        warn=mocker.Mock(), info=mocker.Mock(), debug=mocker.Mock(), error=mocker.Mock())


@pytest.fixture
def get_scorer_kw():
    return dict(runs=10, rands=10, limit=100)


def test_init(lextstat_factory, test_data, mocker, log, tmppath):
    lextstat_factory({0: ['ID', 'doculect', 'concept', 'IPA'],
                    1: ['1', 'deu', 'hand', 'hant']}, model='sca')
    ls = lextstat_factory({0: ['ID', 'doculect', 'concept', 'IPA'],
                         1: ['1', 'deu', 'hand', 'hant']})
    assert 'lexstat' in repr(ls)
    lextstat_factory(ls)
    lextstat_factory({0: ['ID', 'doculect', 'concept', 'tokens'],
                    1: ['1', 'deu', 'hand', ['h', 'a', 'n', 't']]})
    with pytest.raises(AssertionError):
        LexStat({0: ['ID', 'doculect', 'concept'], 1: ['1', 'deu', 'hand']})
    lextstat_factory(str(test_data / 'phybo.qlc'), check=True)
    mocker.patch('lingpy.compare.lexstat.log', log)
    lextstat_factory(str(test_data / 'KSL.qlc'), check=True)
    assert log.info.called
    error_log = tmppath / 'errors'
    mocker.patch('lingpy.util.confirm', mocker.Mock(return_value=True))
    lex = lextstat_factory(
        {
            0: ['ID', 'doculect', 'concept', 'IPA', 'tokens'],
            1: ['1', 'deu', 'hand', 'hand', ['']],
            2: ['2', 'eng', 'hand', 'hand', ['abc']],
            3: ['3', 'xyz', 'hand', 'hund', 'h u n d'],
        },
        check=True, errors='%s' % error_log)
    assert error_log.exists()
    assert lex.filename.endswith('_cleaned.tsv')
    assert pathlib.Path(lex.filename).exists()
    pathlib.Path(lex.filename).unlink()
    assert len(lex._meta['errors']) == 2


def test_init2(lex, test_data):
    freqs = lex.freqs['Hawaiian']
    seq = {'5.W.C': 19, '5.I.V': 87, '5.Y.V': 75, '5.U.V': 87}

    for char, n in seq.items():
        assert freqs[char] == n

    assert len(lex.chars) == 187
    assert len(lex.rchars) == 35

    for name in 'bscorer rscorer pairs'.split():
        obj = jsonlib.load(test_data / 'KSL.{}.json'.format(name))
        if name != 'pairs':
            assert getattr(lex, name).matrix == obj
        else:
            for key, values in lex.pairs.items():
                values = set(values)
                ovalues = set(tuple(v) for v in obj['---'.join(key)])
                if key == 'pairs':
                    assert values == ovalues  # FIXME: This is not hit in tests! Why?


def test_init3(test_data, lextstat_factory):  # with kw check=True
    bad_file = test_data / 'bad_file.tsv'
    with pytest.raises(ValueError):
        LexStat(str(bad_file))
    ls = lextstat_factory(str(bad_file), check=True, apply_checks=True)
    assert hasattr(ls, 'errors')
    cleaned = bad_file.parent.joinpath(bad_file.name + '_cleaned.tsv')
    assert cleaned.exists()
    cleaned.unlink()
    with pytest.raises(ValueError):
        LexStat({0: ['concept', 'language', 'ipa']})


def test_getitem(lex):
    assert lex['xyz'] is None


def test_get_scorer(lex, mocker, get_scorer_kw, log):
    lex.get_scorer(**get_scorer_kw)
    assert hasattr(lex, "cscorer")
    mocker.patch('lingpy.compare.lexstat.log', log)
    lex.get_scorer(**get_scorer_kw)
    assert log.warning.called
    del lex.cscorer
    lex.get_scorer(**get_scorer_kw)
    lex.get_scorer(method='markov', **get_scorer_kw)


def test_cluster(lex, mocker, get_scorer_kw):
    lex.get_scorer(**get_scorer_kw)
    lex.cluster(method="lexstat", threshold=0.7)
    lex.cluster(method="edit-dist", threshold=0.7)
    lex.cluster(method="turchin", threshold=0.7)
    with pytest.raises(ValueError):
        lex.cluster(method="fuzzy")
    mocker.patch('lingpy.basic.parser.confirm', mocker.Mock(return_value=True))
    lex.cluster(method="sca", guess_threshold=True, gt_mode='nulld')
    assert all(x in lex.header for x in 'scaid lexstatid editid turchinid'.split())


def test_align_pairs(lex):
    assert not lex.align_pairs('English', 'German', method='sca', pprint=False)
    assert lex.align_pairs(1, 2, method='sca', pprint=False)[-1] > 0.5


def test__get_matrices(lex):
    matrix = list(lex._get_matrices(concept="hand", method="sca"))[0]
    assert len(matrix) == 7

    matrix = list(lex._get_matrices(concept="hand", method="turchin"))[0]
    assert matrix[0][1] == 1


def test_get_subset(test_data, lex):
    lex.get_subset([])
    assert [v for v in lex.subsets.values() if v] == []
    pairs = jsonlib.load(test_data / 'KSL.pairs.json')
    assert sorted('---'.join(k) for k in lex.subsets.keys()) ==\
        sorted(pairs.keys())


def test_get_distances(lex, get_scorer_kw):
    lex.get_scorer(**get_scorer_kw)
    lex.get_random_distances()
    lex.get_distances()
    lex.get_distances(method='edit-dist')
    lex.get_distances(aggregate=False)


def test_get_frequencies(lex):
    f = lex.get_frequencies('sounds')
    assert len(f) == lex.width

    f = lex.get_frequencies('sounds', aggregated=True)
    tokens = []
    for k in lex:
        for t in lex[k, 'tokens']:
            tokens += [t]
    assert len(f) == len(set(tokens))

    d = lex.get_frequencies('diversity', ref='cogid')
    assert isinstance(d, float)

    w = lex.get_frequencies('wordlength')
    assert len(w) == lex.width

    w = lex.get_frequencies('wordlength', aggregated=True)
    assert isinstance(w, float)


def test_output(lex, tmppath):
    lex.output('csv', filename=str(tmppath /'test_lexstat'))
    lex.output('scorer', filename=str(tmppath / 'test_lexstat'))


def test_correctness(lextstat_factory):
    lex = lextstat_factory({
        0: ['ID', 'doculect', 'concept', 'IPA'],
        1: ['1', 'deu', 'hand', 'hand'],
        2: ['2', 'eng', 'hand', 'hand'],
        3: ['3', 'xyz', 'hand', 'xyz']})
    lex.cluster(ref='cogid', method='sca', threshold=0.5)
    assert lex[1, 'cogid'] == lex[2, 'cogid']

    rc(schema='asjp')
    lex = lextstat_factory({
        0: ['ID', 'concept', 'ipa', 'doculect'],
        1: ['5424', 'Abend::N', 'swar', 'FRA'],
        2: ['5425', 'Abend::N', 'sware', 'FRA'],
        3: ['5426', 'Abend::N', 'sear3', 'RON'],
        4: ['5427', 'Abend::N', 'ivniN', 'ENG'],
        5: ['5428', 'Abend::N', 'noyt3', 'POR'],
        6: ['5429', 'Abend::N', 'tardi5a', 'POR'],
        7: ['5430', 'Abend::N', 'afd3n', 'DAN'],
    })
    lex.cluster(method='sca', threshold=0.5, ref='cogid')
    assert lex[1, 'cogid'], lex[2, 'cogid'] == lex[3, 'cogid']
    rc(schema='ipa')
