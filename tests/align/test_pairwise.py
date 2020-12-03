import pytest

from lingpy.align import (
    Pairwise, pw_align, nw_align, sw_align, we_align, structalign, turchin,
    edit_dist
)
from lingpy.data.model import Model


@pytest.fixture
def pair():
    return Pairwise('waldemar', 'vladimir')


def test_basics(pair):
    assert 'waldemar' in str(pair)
    assert 'waldemar' in repr(pair)
    assert len(pair) == 1
    pair.align()
    assert '-' in str(pair)
    for d in 'wcta':
        assert len(pair[0, d][0]) >= len('waldemar')


def test_align(pair):
    for mode in ['global', 'local', 'overlap', 'dialign']:
        pair.align(mode=mode, distance=True, pprint=True)
        assert '-' in ''.join(pair.alignments[0][1])


def test_pw_align():
    for mode in ['global', 'local', 'overlap', 'dialign']:
        alms1 = pw_align('waldemar', 'vladimir', mode=mode, distance=True)
        alms2 = pw_align('waldemar', 'vladimir', mode=mode, scale=0.5,
                         distance=False)

        assert alms1[-1] != alms2[-1]

    assert pw_align('waldemar', 'vladimir', mode='local')[2] == 1

    scorer = {('a', 'a'): 1}
    assert pw_align('a', 'a', scorer=scorer)[2] == 1


def test_nw_align():
    assert nw_align('waldemar', 'vladimir')[2] == -1
    with pytest.raises(ValueError):
        nw_align(1, 2)


def test_sw_align():
    assert sw_align('waldemar', 'wald')[2] == 4
    with pytest.raises(ValueError):
        sw_align(1, 2)


def test_we_align():
    assert we_align('waldemar', 'marwalde')[1][2] == 3
    with pytest.raises(ValueError):
        we_align(1, 2)


def test_structalign():
    assert len(structalign('waldemar', 'vladimir')[0]) == 3


def test_turchin():
    assert turchin('waldemar', 'vladimir') == 0
    assert turchin('waldemar', 'vladimir', model=Model('dolgo')) == 0


def test_editdist():
    assert edit_dist('waldemar', 'vladimir', restriction="cv") == 5
