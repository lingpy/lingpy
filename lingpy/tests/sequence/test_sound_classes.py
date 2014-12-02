# *-* coding: utf-8 *-* 
# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-12-02 15:49
# modified : 2014-12-02 15:49
"""
Test sound class functions.
"""

__author__="Johann-Mattis List"
__date__="2014-12-02"

# modified : 2013-11-12 12:53
"""
Test the SCA module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-12"

import os
import unittest
from lingpy.sequence.sound_classes import ipa2tokens, token2class, \
        tokens2class, prosodic_string, prosodic_weights, class2tokens, pid,\
        check_tokens, get_all_ngrams, sampa2uni, bigrams, trigrams, fourgrams,\
        get_n_ngrams, pgrams
from lingpy import rc


def test_ipa2tokens():

    seq = 'tʰɔxtər'

    assert len(ipa2tokens(seq)) != len(list(seq))

def test_token2class():

    seq = 'tʰ ɔ x t ə r'.split(' ')

    assert token2class(seq[0], rc('dolgo')) == 'T'

def test_tokens2class():

    seq = 'tʰ ɔ x t ə r'.split(' ')

    assert tokens2class(seq, rc('dolgo')) == list('TVKTVR')

def test_prosodic_string():

    seq = 'tʰ ɔ x t ə r'.split(' ')

    assert prosodic_string(seq) == 'AXMBYN'

def test_prosodic_weights():

    seq = 'tʰ ɔ x t ə r'.split(' ')

    assert prosodic_weights(prosodic_string(seq))[0] == 2
    assert prosodic_weights(prosodic_string(seq))[-1] == 0.8

def test_class2tokens():

    classes = 'T-VKTV-R'
    tokens = 'tʰ ɔ x t ə r'.split(' ')

    out = class2tokens(classes, tokens)
    assert out[1] == '-' and out[-2] == '-'

def test_pid():

    assert pid('mattis', 'maTTIs') == 0.5

def test_check_tokens():

    assert check_tokens('th o x T e r'.split(' '))[0] == (3,'T')

def test_get_all_ngrams():

    assert get_all_ngrams('ab') == ['ab', 'a', 'b']

def test_sampa2uni():

    seq = 'tʰɔxtər'
    
    assert sampa2uni('t_hOxt@r') == seq

def bigrams():

    assert bigrams('ab')[0] == ('#','a') and bigrams('ab')[-1] == ('b','$')

def trigrams():

    assert bigrams('ab')[0] == ('#', '#' ,'a') and bigrams('ab')[-1] == ('b','$', '$')

def fourgrams():

    assert fourgrams('ab')[0] == ('#','#','#','a')
    assert fourgrams('ab')[-1] == ('a','$','$','$')

def test_pgrams():

    assert pgrams('ab')[0] == ('a','X')
    assert pgrams('ab')[-1] == ('b','N')
