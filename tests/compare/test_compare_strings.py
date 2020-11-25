import pytest

from lingpy.compare import strings


@pytest.fixture
def a():
    return 'levenshtein'


@pytest.fixture
def b():
    return 'lveenshtiene'


def test_ldn(a, b):
    _ = strings.ldn(a, b, normalized=True)
    d = strings.ldn(a, b, normalized=False)
    assert int(d) == 5


def test_ldn_swap(a, b):
    _ = strings.ldn_swap(a, b, normalized=True)
    d = strings.ldn_swap(a, b, normalized=False)
    assert int(d) == 3


def test_bidist1(a, b):
    _ = strings.bidist1(a, b, normalized=True)
    d = strings.bidist1(a, b, normalized=False)
    assert int(d) == 8


def test_tridist1(a, b):
    _ = strings.tridist1(a, b, normalized=True)
    d = strings.tridist1(a, b, normalized=False)
    assert int(d) == 10


def test_bidist2(a, b):
    _ = strings.bidist2(a, b, normalized=True)
    d = strings.bidist2(a, b, normalized=False)
    assert int(d) == 3


def test_tridist2(a, b):
    _ = strings.tridist2(a, b, normalized=True)
    d = strings.tridist2(a, b, normalized=False)
    assert int(d) == 2


def test_bidist3(a, b):
    _ = strings.bidist3(a, b, normalized=True)
    d = strings.bidist3(a, b, normalized=False)
    assert int(d) == 5


def test_tridist3(a, b):
    _ = strings.tridist3(a, b, normalized=True)
    d = strings.tridist3(a, b, normalized=False)
    assert int(d) == 5


def test_dice(a, b):
    _ = strings.dice(a, b, normalized=True)
    d = strings.dice(a, b, normalized=False)
    assert int(d) == 11


def test_lcs(a, b):
    _ = strings.lcs(a, b, normalized=True)
    d = strings.lcs(a, b, normalized=False)
    assert int(d) == 3


def test_bisim1(a, b):
    _ = strings.bisim1(a, b, normalized=True)
    d = strings.bisim1(a, b, normalized=False)
    assert int(d) == 6


def test_trisim1(a, b):
    _ = strings.trisim1(a, b, normalized=True)
    d = strings.trisim1(a, b, normalized=False)
    assert int(d) == 8


def test_bisim2(a, b):
    _ = strings.bisim2(a, b, normalized=True)
    d = strings.bisim2(a, b, normalized=False)
    assert int(d) == 3


def test_trisim2(a, b):
    _ = strings.trisim2(a, b, normalized=True)
    d = strings.trisim2(a, b, normalized=False)
    assert int(d) == 2


def test_bisim3(a, b):
    _ = strings.bisim3(a, b, normalized=True)
    d = strings.bisim3(a, b, normalized=False)
    assert int(d) == 4


def test_trisim3(a, b):
    _ = strings.trisim3(a, b, normalized=True)
    d = strings.trisim3(a, b, normalized=False)
    assert int(d) == 4


def test_jcd(a, b):
    _ = strings.jcd(a, b, normalized=True)
    d = strings.jcd(a, b, normalized=False)
    assert int(d) == 11


def test_jcdn(a, b):
    _ = strings.jcdn(a, b, normalized=True)
    d = strings.jcdn(a, b, normalized=False)
    assert int(d) == 15


def test_prefix(a, b):
    _ = strings.prefix(a, b, normalized=True)
    d = strings.prefix(a, b, normalized=False)
    assert int(d) == 5


def test_xdice(a, b):
    _ = strings.xdice(a, b, normalized=True)
    d = strings.xdice(a, b, normalized=False)
    assert int(d) == 9


def test_trigram(a, b):
    _ = strings.trigram(a, b, normalized=True)
    d = strings.trigram(a, b, normalized=False)
    assert int(d) == 13


def test_ident(a, b):
    d = strings.ident(a, a)
    assert int(d) == 1


def test_xxdice(a, b):
    _ = strings.xxdice(a, b, normalized=True)
    _ = strings.xxdice('aa', '')
    d = strings.xxdice(a, b, normalized=False)
    assert int(d) == 11
