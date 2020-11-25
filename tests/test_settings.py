import pytest

from lingpy.settings import rc


def test_rc_get():
    assert rc(rval='x', rcParams_={'x': 5}) == 5


@pytest.mark.parametrize(
    'key,val,test',
    [
        ('schema', 'ipa', lambda d: 'dolgo' in d),
        ('schema', 'el', lambda d: not d['merge_vowels']),
        ('custom', 5, lambda d: d['custom'] == 5),
        ('fn', 'x', lambda d: 'filename' in d),
    ]
)
def test_rc(key, val, test):
    d = {}
    rc(rcParams_=d, **{key: val})
    assert test(d)
