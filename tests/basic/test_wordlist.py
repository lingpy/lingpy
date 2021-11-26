"""
Test wordlist module.
"""
import pytest

from lingpy import Wordlist


@pytest.fixture
def wordlist(test_data):
    return Wordlist(str(test_data / 'KSL.qlc'))


@pytest.fixture
def wordlist2(test_data):
    return Wordlist(str(test_data / 'good_file.tsv'))


def test___len__(wordlist):
    assert len(wordlist) == 1400


def test_calculate(wordlist):
    wordlist.calculate('dst')

    assert hasattr(wordlist, 'distances')
    assert sum([wordlist.distances[x][x] for x in
                range(wordlist.width)]) == 0

    wordlist.calculate('tree')
    assert str(wordlist.tree).endswith(';')

    assert sorted(wordlist.tree.taxa) == sorted(wordlist.cols)

    wordlist.calculate('groups')
    assert hasattr(wordlist, 'groups')
    assert type(wordlist.groups) == dict


def test_coverage(wordlist):
    wordlist.coverage()
    wordlist.coverage(stats='ratio')
    wordlist.coverage(stats='mean')


def test_get_list(wordlist, wordlist2):
    ger_l = wordlist.get_list(doculect='German', entry='ipa', flat=True)
    ger_d = wordlist.get_dict(col='German', entry='ipa')
    ger_t = wordlist.get_list(doculect='German', entry="ipa")

    assert sorted(ger_l) == sorted([v[0] for v in ger_d.values()])
    assert sorted(ger_t) == sorted(ger_l)

    hand1 = wordlist.get_list(concept="hand", entry="ipa", flat=True)
    hand2 = wordlist.get_dict(row="hand", entry="ipa")
    assert sorted(hand1) == sorted([v[0] for v in hand2.values()])

    # test for synonym lines, which are flattened
    assert wordlist2.get_list(concept='hand', entry="language",
                                   flat=True).count('l6') == 2
    nonflat = wordlist2.get_list(concept="hand", entry="language")
    assert nonflat[0][-1] == nonflat[1][-1]
    assert len(wordlist2.get_list(col="l1", entry="concept")) == 3
    assert len(wordlist2.get_list(col="l1", flat=True,
                                       entry="concept")) == 2

    with pytest.raises(ValueError):
        wordlist2.get_list(col="l1", row="hand")
    with pytest.raises(ValueError):
        wordlist2.get_list()
    with pytest.raises(ValueError):
        wordlist.get_list(**{"row": "Hand"})


def test_get_dict(wordlist):
    ger_d = wordlist.get_dict(col='German')

    assert sorted(ger_d.keys()) == sorted(wordlist.rows)
    with pytest.raises(ValueError):
        wordlist.get_dict(**{"row": "Hand"})


def test_renumber(wordlist):
    wordlist.renumber('cogid', 'dummy')

    ger1 = wordlist.get_list(col='German', entry='cogid', flat=True)
    ger2 = wordlist.get_list(col='German', entry='dummy', flat=True)

    assert len(set(ger1)) == len(set(ger2))
    assert sum([1 for x in ger2 if type(x) == int]) == len(ger2)


def test_get_entries(wordlist):
    ger = wordlist.get_entries('cogid')

    assert len(ger) == wordlist.height
    assert len(ger[0]) == wordlist.width


def test_get_etymdict(wordlist):
    etd1 = wordlist.get_etymdict(ref='cogid', entry='ipa',
                                      modify_ref=False)
    etd2 = wordlist.get_etymdict(ref='cogid', entry='ipa',
                                      modify_ref=abs)


    assert (len(etd1) > len(etd2) and
            len(set([abs(x) for x in etd1])) == len(etd2))
    assert len([x for x in etd2 if x < 0]) == 0

    # test iter_cognates
    for cogid, idxs in wordlist.iter_cognates("cogid"):
        assert len(idxs) == 1

    # make "fuzzy" cognate sets
    wordlist.add_entries('fuzzyid', 'cogid', lambda x: [x])

    etd3 = wordlist.get_etymdict(ref='fuzzyid', entry='ipa', modify_ref=False)
    etd4 = wordlist.get_etymdict(ref='fuzzyid', entry='ipa', modify_ref=abs)
    for key in etd1:
        assert etd1[key] == etd3[key]
    for key in etd2:
        assert etd2[key] == etd4[key]


def test_get_paps(wordlist):
    paps = wordlist.get_paps(ref="cogid", modify_ref=abs)
    cogs = wordlist.get_etymdict(ref="cogid", modify_ref=abs)

    for key in cogs:
        assert abs(key) in paps


def test_output(tmp_path, wordlist):
    fn = str(tmp_path / 'test')
    for fmt in 'tsv taxa tre dst starling paps.nex paps.csv' \
               'separated multistate.nex groups'.split():
        kw = {'ref': 'word'} if fmt == 'starling' else {}
        wordlist.output(fmt, filename=fn, **kw)

        if fmt == 'starling':
            wordlist.output(fmt, filename=fn, cognates='cogid', **kw)
        if fmt == 'tsv':
            kw['subset'] = True
            wordlist.output(fmt, filename=fn, cols=[], rows={}, **kw)
            wordlist.output(fmt, filename=fn,
                                 cols=sorted(wordlist.header)[:2],
                                 rows=dict(ID=" > 10"), **kw)


def test_export(tmp_path, wordlist):
    fn = str(tmp_path / 'test')
    for fmt in 'txt tex html'.split():
        wordlist.export(fmt, filename=fn)


def test_get_wordlist(test_data):
    from lingpy.basic.wordlist import get_wordlist
    wl1 = get_wordlist(str(test_data / 'mycsvwordlist.csv'))
    wl2 = get_wordlist(str(test_data / 'mycsvwordlistwithoutids.csv'))
    assert wl1.height == wl2.height
    for k in wl1:
        assert wl1[k, 'concept'] == wl2[k, 'concept']
