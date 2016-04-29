from __future__ import print_function, division, unicode_literals
import os
from pathlib import Path

from mock import patch, Mock
from nose.tools import assert_raises
from lingpy import LexStat, rc
from lingpy.compare.lexstat import char_from_charstring, get_score_dict
from lingpy.util import jsonload
from lingpy.tests.util import test_data, WithTempDir, get_log

def test_char_from_charstring():
    assert char_from_charstring('a.b.c') == "b"
    assert char_from_charstring('a.b') == "a"
    assert_raises(ValueError, char_from_charstring, "a")

def test_get_score_dict():
    chars = ["1.A.-", "2.B.-"]
    model = rc("sca")
    sd = get_score_dict(chars, model)
    assert sd['A','B'] == -22.5
    

class TestLexStat(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.lex = self._make_one(test_data('KSL.qlc'))
        self.log = get_log()
        self.get_scorer_kw = dict(runs=10, rands=10, limit=100)

    def _make_one(self, *args, **kw):
        kw.setdefault('errors', self.tmp_path('errors.log').as_posix())
        return LexStat(*args, **kw)

    def test_init(self):
        self._make_one({0: ['ID', 'doculect', 'concept', 'IPA'],
                        1: ['1', 'deu', 'hand', 'hant']}, model='sca')
        ls = self._make_one({0: ['ID', 'doculect', 'concept', 'IPA'],
                             1: ['1', 'deu', 'hand', 'hant']})
        self.assertIn('lexstat', repr(ls))
        self._make_one(ls)
        self._make_one({0: ['ID', 'doculect', 'concept', 'tokens'],
                        1: ['1', 'deu', 'hand', 'hant']})
        self.assertRaises(AssertionError, LexStat, {0: ['ID', 'doculect',
                                                        'concept'],
                                                    1: ['1', 'deu', 'hand']})
        self._make_one(test_data('phybo.qlc'), check=True)
        with patch('lingpy.compare.lexstat.log', self.log):
            self._make_one(test_data('KSL.qlc'), check=True)
            assert self.log.info.called
        error_log = self.tmp_path('errors')
        with patch('lingpy.util.confirm', Mock(return_value=True)):
            lex = self._make_one({
                0: ['ID', 'doculect', 'concept', 'IPA', 'tokens'],
                1: ['1', 'deu', 'hand', 'hand', ['']],
                2: ['2', 'eng', 'hand', 'hand', ['abc']],
                3: ['3', 'xyz', 'hand', 'hund', 'h u n d'],
            }, check=True, errors='%s' % error_log)
            assert error_log.exists()
            self.assertTrue(lex.filename.endswith('_cleaned.tsv'))
            self.assertTrue(os.path.exists(lex.filename))
            os.remove(lex.filename)
            self.assertEquals(len(lex._meta['errors']), 2)

    def test_init2(self):
        freqs = self.lex.freqs['Hawaiian']
        for char, n in {'5.W.C': 19, '5.I.V': 87, '5.Y.V': 75, '5.U.V': 87}.items():
            self.assertEquals(freqs[char], n)
        self.assertEquals(len(self.lex.chars), 187)
        self.assertEquals(len(self.lex.rchars), 35)

        self.maxDiff = None

        for name in 'bscorer rscorer pairs'.split():
            obj = jsonload(test_data('KSL.%s.json' % name))
            if name != 'pairs':
                self.assertEquals(getattr(self.lex, name).matrix, obj)
            else:
                for key, values in self.lex.pairs.items():
                    values = set(values)
                    ovalues = set(tuple(v) for v in obj['---'.join(key)])
                    if name != 'pairs':
                        self.assertEquals(values, ovalues)

    def test_init3(self):  # with kw check=True
        bad_file = Path(test_data('bad_file.tsv'))
        assert_raises(ValueError, LexStat, bad_file.as_posix())
        ls = self._make_one(bad_file.as_posix(), check=True, apply_checks=True)
        assert hasattr(ls, 'errors')
        cleaned = bad_file.parent.joinpath(bad_file.name + '_cleaned.tsv')
        self.assertTrue(cleaned.exists())
        os.remove(cleaned.as_posix())
        assert_raises(ValueError, LexStat, {0: ['concept', 'language', 'ipa']})

    def test_getitem(self):
        self.assertIsNone(self.lex['xyz'])

    def test_get_scorer(self):
        self.lex.get_scorer(**self.get_scorer_kw)
        assert hasattr(self.lex, "cscorer")
        with patch('lingpy.compare.lexstat.log', self.log):
            self.lex.get_scorer(**self.get_scorer_kw)
            assert self.log.warn.called
        del self.lex.cscorer
        self.lex.get_scorer(**self.get_scorer_kw)
        self.lex.get_scorer(method='markov', **self.get_scorer_kw)

    def test_cluster(self):
        self.lex.get_scorer(**self.get_scorer_kw)
        self.lex.cluster(method="lexstat", threshold=0.7)
        self.lex.cluster(method="edit-dist", threshold=0.7)
        self.lex.cluster(method="turchin", threshold=0.7)
        self.assertRaises(ValueError, self.lex.cluster, method="fuzzy")
        with patch('lingpy.basic.parser.confirm', Mock(return_value=True)):
            self.lex.cluster(method="sca", guess_threshold=True, gt_mode='nulld')

        assert 'scaid' in self.lex.header \
               and 'lexstatid' in self.lex.header \
               and 'editid' in self.lex.header \
               and 'turchinid' in self.lex.header

    def test_align_pairs(self):
        assert not self.lex.align_pairs('English', 'German', method='sca', pprint=False) 
        assert self.lex.align_pairs(1, 2, method='sca', pprint=False)[-1] > 0.5
    
    def test__get_matrices(self):

        matrix = list(self.lex._get_matrices(concept="hand", method="sca"))[0]
        assert len(matrix) == 7
        
        matrix = list(self.lex._get_matrices(concept="hand",
            method="turchin"))[0]
        assert matrix[0][1] == 1

    def test_get_subset(self):
        self.lex.get_subset([])
        self.assertEquals([v for v in self.lex.subsets.values() if v], [])
        pairs = jsonload(test_data('KSL.pairs.json'))
        self.assertEquals(
            sorted('---'.join(k) for k in self.lex.subsets.keys()),
            sorted(pairs.keys()))

    def test_get_distances(self):
        self.lex.get_scorer(**self.get_scorer_kw)
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

    def test_correctness(self):
        lex = self._make_one({
            0: ['ID', 'doculect', 'concept', 'IPA'],
            1: ['1', 'deu', 'hand', 'hand'],
            2: ['2', 'eng', 'hand', 'hand'],
            3: ['3', 'xyz', 'hand', 'xyz']})
        lex.cluster(ref='cogid', method='sca', threshold=0.5)
        self.assertEquals(lex[1, 'cogid'], lex[2, 'cogid'])

        rc(schema='asjp')
        lex = self._make_one({
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
        self.assertEquals(lex[1, 'cogid'], lex[2, 'cogid'], lex[3, 'cogid'])
        rc(schema='ipa')
