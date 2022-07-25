"""
Test wordlist module.
"""
import pytest

from lingpy import Wordlist, Alignments
from lingpy.basic.ops import wl2dst, wl2qlc, tsv2triple, triple2tsv, \
    calculate_data, wl2multistate, coverage, clean_taxnames, wl2dict


@pytest.fixture
def wordlist(test_data):
    return Wordlist(str(test_data / 'KSL.qlc'))


@pytest.fixture
def wordlist2(test_data):
    return Wordlist(str(test_data / 'KSL4.qlc'))


def test_iter_rows(wordlist):
    from lingpy.basic.ops import iter_rows
    assert len(list(
        iter_rows(wordlist, 'concept', 'doculect'))[0]) == 3


def test_wl2dict(wordlist):
    _ = wl2dict(wordlist, dict(s1=['concept', '{0}'], s2=['cogid', '{0}']), [('ipa', '{0}')])


def test_wl2dst(wordlist):
    res = wl2dst(wordlist, mode='jaccard')
    assert isinstance(res, list)
    res = wl2dst(wordlist, mode='jaccard', refB='glossid')
    assert isinstance(res, list)

    _ = wl2dst(wordlist, mode='swadesh')
    _ = wl2dst(wordlist, mode='shared')
    _ = wl2dst(wordlist, mode='swadesh', ignore_missing=True)

    # trigger zero-division-warning in wl2dst
    tmp = Wordlist({
        0: ['doculect', 'concept', 'counterpart', 'cogid'],
        1: ['l1', 'hand', 'hand', '1'],
        2: ['l2 - a (taxon) name)', 'hand', 'hand', '2'],
        3: ['l3', 'foot', 'foot', '3']
    })
    dst = wl2dst(tmp)
    assert dst[0][2] == 1


def test_wl2qlc(tmp_path, test_data, wordlist):
    stamp = 'test-stamp'
    out = tmp_path / 'test'

    wl2qlc(wordlist.header, wordlist._data, filename=str(out), stamp=stamp)
    out = tmp_path / 'test.qlc'
    assert out.read_text(encoding='utf8').endswith(stamp)
    
    # use pathlib instance
    wl2qlc(wordlist.header, wordlist._data, filename=out, stamp=stamp)
    out = tmp_path / 'test.qlc'
    assert out.read_text(encoding='utf8').endswith(stamp)

    # load a wordlist with alignments and otuput it as string with msapairs
    tmp = Alignments(str(test_data / 'good_file.tsv'), ref='cogid')
    tmp.align(ref="cogid")

    wl2qlc(tmp.header, tmp._data, meta=tmp._meta, filename=str(out), stamp='stampo', ignore=[])
    tmp.get_consensus(ref="cogid")

    wl2qlc([h.upper()
            for h in sorted(tmp.header, key=lambda x: tmp.header[x])],
           tmp._data, meta=tmp._meta, filename=out.as_posix(),
           stamp='stampo', ignore=[], formatter="doculect,concept")
    wl2qlc([h.upper()
            for h in sorted(tmp.header, key=lambda x: tmp.header[x])],
           tmp._data, meta=tmp._meta, filename=out.as_posix(),
           stamp='stampo', ignore=[], formatter="doculect")


def test_tsv2triple(tmp_path, wordlist):
    out = tmp_path / 'test'
    triples = tsv2triple(wordlist, str(out))
    assert isinstance(triple2tsv(str(out)), list)
    assert isinstance(triple2tsv(triples, output='dict'), dict)


def test_calculate_data(wordlist):
    for data in ['tree', 'dst', 'cluster']:
        calculate_data(wordlist, data)
        calculate_data(wordlist, data, mode='shared')


def test_wl2multistate(wordlist2):
    res = wl2multistate(wordlist2, 'cogid', '?')
    # the result must be a matrix.
    assert isinstance(res, list)
    assert len(set(len(row) for row in res)) == 1


def test_coverage(wordlist):
    res = coverage(wordlist)
    assert res['Turkish'] == 200


def test_clean_taxnames():
    tmp = Wordlist({
        0: ['doculect', 'concept', 'counterpart'],
        1: ['l1', 'hand', 'hand'],
        2: ['l2 - a (taxon) name)', 'hand', 'hand']
    })

    clean_taxnames(tmp)
    assert tmp.cols[-1] == 'l2___a_taxon_name'


def test_renumber(test_data):
    tmp = Wordlist(str(test_data / 'good_file.tsv'))
    tmp.renumber('cogid', 'newcogid')
    assert 'newcogid' in tmp.header
    tmp.renumber('mock')
    assert 'mockid' in tmp.header
