from unittest import TestCase

from lingpy.compare import strings


class TestStrings(TestCase):
    def setUp(self):
        self.a = 'levenshtein'
        self.b = 'lveenshtiene'

    def test_ldn(self):
        _ = strings.ldn(self.a, self.b, normalized=True)
        d = strings.ldn(self.a, self.b, normalized=False)
        assert int(d) == 5

    def test_ldn_swap(self):
        _ = strings.ldn_swap(self.a, self.b, normalized=True)
        d = strings.ldn_swap(self.a, self.b, normalized=False)
        assert int(d) == 3

    def test_bidist1(self):
        _ = strings.bidist1(self.a, self.b, normalized=True)
        d = strings.bidist1(self.a, self.b, normalized=False)
        assert int(d) == 8

    def test_tridist1(self):
        _ = strings.tridist1(self.a, self.b, normalized=True)
        d = strings.tridist1(self.a, self.b, normalized=False)
        assert int(d) == 10

    def test_bidist2(self):
        _ = strings.bidist2(self.a, self.b, normalized=True)
        d = strings.bidist2(self.a, self.b, normalized=False)
        assert int(d) == 3

    def test_tridist2(self):
        _ = strings.tridist2(self.a, self.b, normalized=True)
        d = strings.tridist2(self.a, self.b, normalized=False)
        assert int(d) == 2

    def test_bidist3(self):
        _ = strings.bidist3(self.a, self.b, normalized=True)
        d = strings.bidist3(self.a, self.b, normalized=False)
        assert int(d) == 5

    def test_tridist3(self):
        _ = strings.tridist3(self.a, self.b, normalized=True)
        d = strings.tridist3(self.a, self.b, normalized=False)
        assert int(d) == 5

    def test_dice(self):
        _ = strings.dice(self.a, self.b, normalized=True)
        d = strings.dice(self.a, self.b, normalized=False)
        assert int(d) == 11

    def test_lcs(self):
        _ = strings.lcs(self.a, self.b, normalized=True)
        d = strings.lcs(self.a, self.b, normalized=False)
        assert int(d) == 3

    def test_bisim1(self):
        _ = strings.bisim1(self.a, self.b, normalized=True)
        d = strings.bisim1(self.a, self.b, normalized=False)
        assert int(d) == 6

    def test_trisim1(self):
        _ = strings.trisim1(self.a, self.b, normalized=True)
        d = strings.trisim1(self.a, self.b, normalized=False)
        assert int(d) == 8

    def test_bisim2(self):
        _ = strings.bisim2(self.a, self.b, normalized=True)
        d = strings.bisim2(self.a, self.b, normalized=False)
        assert int(d) == 3

    def test_trisim2(self):
        _ = strings.trisim2(self.a, self.b, normalized=True)
        d = strings.trisim2(self.a, self.b, normalized=False)
        assert int(d) == 2

    def test_bisim3(self):
        _ = strings.bisim3(self.a, self.b, normalized=True)
        d = strings.bisim3(self.a, self.b, normalized=False)
        assert int(d) == 4

    def test_trisim3(self):
        _ = strings.trisim3(self.a, self.b, normalized=True)
        d = strings.trisim3(self.a, self.b, normalized=False)
        assert int(d) == 4

    def test_jcd(self):
        _ = strings.jcd(self.a, self.b, normalized=True)
        d = strings.jcd(self.a, self.b, normalized=False)
        assert int(d) == 11

    def test_jcdn(self):
        _ = strings.jcdn(self.a, self.b, normalized=True)
        d = strings.jcdn(self.a, self.b, normalized=False)
        assert int(d) == 15

    def test_prefix(self):
        _ = strings.prefix(self.a, self.b, normalized=True)
        d = strings.prefix(self.a, self.b, normalized=False)
        assert int(d) == 5

    def test_xdice(self):
        _ = strings.xdice(self.a, self.b, normalized=True)
        d = strings.xdice(self.a, self.b, normalized=False)
        assert int(d) == 9

    def test_trigram(self):
        _ = strings.trigram(self.a, self.b, normalized=True)
        d = strings.trigram(self.a, self.b, normalized=False)
        assert int(d) == 13

    def test_ident(self):
        d = strings.ident(self.a, self.a)
        assert int(d) == 1

    def test_xxdice(self):
        _ = strings.xxdice(self.a, self.b, normalized=True)
        _ = strings.xxdice('aa', '')
        d = strings.xxdice(self.a, self.b, normalized=False)
        assert int(d) == 11
