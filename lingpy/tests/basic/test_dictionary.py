from __future__ import unicode_literals
# author   : Peter Bouda
# email    : pbouda@cidles.eu
# created  : 2013-08-21 08:43
"""
This contains the test classes and functions for dictionary.py.

"""

__author__="Peter Bouda"
__date__="2013-08-21"

import os
from unittest import TestCase

from lingpy import Dictionary
from lingpy.tests.util import test_data


class TestDictionary(TestCase):

    def setUp(self):
        self.dictionary = Dictionary(test_data('leach1969-67-161.csv'))

    def test___getitem__(self):
        item = self.dictionary[1]
        assert item == ['leach1969/67/1', 'aa', 'Ocaina', 'venir', 'Castellano']

    def test___len__(self):
        assert len(self.dictionary) == 5421

    def test___getattr__(self):
        head_iso = self.dictionary['head_iso']
        assert head_iso == "oca"

    def test_get_tuples(self):
        tuples = self.dictionary.get_tuples()
        assert(tuples[0]) == ('aa', 'venir')

    def test_add_entries(self):
        function = lambda x: x.split(' ')
        self.dictionary.add_entries(
            'tokens',
            'translation',
            function)
        tuples = self.dictionary.get_tuples(['translation', 'tokens'])
        assert tuples[5] == \
            ('los que vinieron (hace tiempo)',
                ['los', 'que', 'vinieron', '(hace', 'tiempo)'])

    def test_tokenize(self):
        self.dictionary.tokenize()
        tuples = self.dictionary.get_tuples(['head', 'tokens'])
        self.assertEqual(tuples[0], ('aa', ['a', 'a']))
