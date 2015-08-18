from __future__ import print_function, division, unicode_literals

from mock import patch, Mock

from lingpy import LexStat
from lingpy.tests.util import test_data, WithTempDir


class TestLexStat(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.lex = LexStat(test_data('KSL.qlc'))

    def test_init(self):
        LexStat(test_data('phybo.qlc'), check=True)
        LexStat(test_data('KSL.qlc'), check=True)

    def test_get_scorer(self):
        self.lex.get_scorer()
        assert hasattr(self.lex, "cscorer")
        del self.lex.cscorer
        self.lex.get_scorer(method='markov')

    def test_cluster(self):
        self.lex.get_scorer()
        self.lex.cluster(method="lexstat", threshold=0.7)
        self.lex.cluster(method="edit-dist", threshold=0.7)
        self.lex.cluster(method="turchin", threshold=0.7)
        self.assertRaises(ValueError, self.lex.cluster, method="fuzzy")
        with patch('lingpy.basic.parser.input', Mock(return_value='y')):
            self.lex.cluster(method="sca", guess_threshold=True, gt_mode='nulld')

        assert 'scaid' in self.lex.header \
            and 'lexstatid' in self.lex.header \
            and 'editid' in self.lex.header \
            and 'turchinid' in self.lex.header

    def test_align_pairs(self):
        self.lex.align_pairs('English', 'German', method='sca')

    def test_get_subset(self):
        self.lex.get_subset([])
        self.assertEquals([v for v in self.lex.subsets.values() if v], [])

    def test_get_distances(self):
        self.lex.get_scorer()
        self.lex.get_random_distances()
        self.lex.get_distances()
        self.lex.get_distances(method='turchin')
        self.lex.get_distances(aggregate=False)

    def test_get_frequencies(self):
        f = self.lex.get_frequencies('sounds')
        assert len(f) == self.lex.width

        f = self.lex.get_frequencies('sounds', aggregated=True)
        tokens = []
        for k in self.lex:
            for t in self.lex[k, 'tokens']:
                tokens += [t]
        assert len(f) == len(set(tokens))

        d = self.lex.get_frequencies('diversity', ref='cogid')
        assert isinstance(d, float)

        w = self.lex.get_frequencies('wordlength')
        assert len(w) == self.lex.width

        w = self.lex.get_frequencies('wordlength', aggregated=True)
        assert isinstance(w, float)

    def test_output(self):
        self.lex.output('csv', filename='%s' % self.tmp_path('test_lexstat'))
        self.lex.output('scorer', filename='%s' % self.tmp_path('test_lexstat'))
