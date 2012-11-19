# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
# Copyright (c) 2011, 2012, Quantitative Language Comparison Team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os, types
import numpy.testing

from qlc.ngram import ngrams_from_graphemes,\
    words_ngrams_matrix_for_graphemes_list

class testNgram(numpy.testing.TestCase):
    """
    The test class for the ngram functions.
    """
    
    def test_ngrams_from_graphemes(self):
        test_word = "ababc"
        test_word_ngrams = ["ab","ba","ab","bc"]
        resulting_ngrams = ngrams_from_graphemes(test_word)
        for ngram in resulting_ngrams:
            assert(ngram in test_word_ngrams)
        
    def test_words_ngrams_matrix_for_graphemes_list(self):
        test_list_words = ["abab","cdab","abcdcd","cd"]
        resulting_matrix = \
            words_ngrams_matrix_for_graphemes_list(test_list_words)
        assert(len(resulting_matrix[0]) == len(test_list_words))
