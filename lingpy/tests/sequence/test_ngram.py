from __future__ import unicode_literals
from unittest import TestCase


class TestNgram(TestCase):
    def test_character_model(self):
        from lingpy.sequence.ngram import character_model

        counts = {
            c: n for c, n, _ in character_model(['abcd', 'abc', 'ab', 'a'], test=True)}
        self.assertEqual(counts, {'a': 4, 'b': 3, 'c': 2, 'd': 1})

    def test_ngrams_from_graphemes(self):
        from lingpy.sequence.ngram import ngrams_from_graphemes

        graphemes = ('#', 'h', '#')
        ngrams_from_graphemes(graphemes, 1)