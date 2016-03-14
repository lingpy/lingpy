from unittest import TestCase

from lingpy.compare import strings


class TestStrings(TestCase):
    def setUp(self):
        self.a = 'levenshtein'
        self.b = 'lveenshtiene'

    def test_ldn(self):
        d = strings.ldn(self.a, self.b)
        assert d == 5

    def test_ldn_swap(self):
        d = strings.ldn_swap(self.a, self.b)
        assert d == 3

    def test_bidist1(self):
        strings.bidist1(self.a, self.b)

    def test_tridist1(self):
        strings.tridist1(self.a, self.b)

    def test_bidist2(self):
        strings.bidist2(self.a, self.b)

    def test_tridist2(self):
        strings.tridist2(self.a, self.b)

    def test_bidist3(self):
        strings.bidist3(self.a, self.b)

    def test_tridist3(self):
        strings.tridist3(self.a, self.b)

    def test_dice(self):
        strings.dice(self.a, self.b)

    def test_lcs(self):
        strings.lcs(self.a, self.b)

    def test_bisim1(self):
        strings.bisim1(self.a, self.b)

    def test_trisim1(self):
        strings.trisim1(self.a, self.b)

    def test_bisim2(self):
        strings.bisim2(self.a, self.b)

    def test_trisim2(self):
        strings.trisim2(self.a, self.b)

    def test_bisim3(self):
        strings.bisim3(self.a, self.b)

    def test_trisim3(self):
        strings.trisim3(self.a, self.b)

    def test_jcd(self):
        strings.jcd(self.a, self.b)

    def test_jcdn(self):
        strings.jcdn(self.a, self.b)

    def test_prefix(self):
        strings.prefix(self.a, self.b)

    def test_xdice(self):
        strings.xdice(self.a, self.b)

    def test_trigram(self):
        strings.trigram(self.a, self.b)

    def test_ident(self):
        strings.ident(self.a, self.b)

    def test_xxdice(self):
        strings.xxdice(self.a, self.b)
