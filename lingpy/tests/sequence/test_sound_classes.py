# *-* coding: utf-8 *-* 
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from ..util import test_data
from nose.tools import assert_raises
from lingpy.sequence.sound_classes import ipa2tokens, token2class, \
        tokens2class, prosodic_string, prosodic_weights, class2tokens, pid,\
        check_tokens, get_all_ngrams, sampa2uni, bigrams, trigrams, fourgrams,\
        get_n_ngrams, pgrams
from lingpy import rc, csv2list

"""
Tests for the SCA module
"""

def test_ipa2tokens():

    seq = 'ˈtʲʰoɔːix_tərp͡f¹¹'

    assert len(ipa2tokens(seq)) != len(list(seq))

    seq = 'ʰto͡i'
    
    assert len(ipa2tokens(seq)) == 2

    seq = 'th o x t a'
    
    assert len(ipa2tokens(seq)) == len(seq.split(' '))

    seq = '# b l a #'
    
    assert len(ipa2tokens(seq)) == len(seq.split(' '))-2

    # now check with all possible data we have so far, but only on cases where
    # tokenization doesn't require the merge_vowels = False flag
    tokens = csv2list(test_data('test_tokenization.tsv'))
    
    for a,b in tokens:
        
        tks = ' '.join(ipa2tokens(a))

        # we check for two variants, since we don't know whether vowels are
        # merged or not in the test data
        assert tks == b

    # now test on smaller set with unmerged vowels 
    tokens = csv2list(test_data('test_tokenization_mv.tsv'))
    
    for a,b in tokens:
        
        tks = ' '.join(ipa2tokens(a, merge_vowels=False, merge_geminates=False))

        # we check for two variants, since we don't know whether vowels are
        # merged or not in the test data
        assert tks == b

    tokens = csv2list(test_data('test_tokenization_nasals.tsv'))
    for a,b in tokens:
        print(tks)
        tks = ' '.join(ipa2tokens(a, merge_vowels=True, merge_geminates=True,
            expand_nasals=True, semi_diacritics='h'))
        assert tks == b

def test_token2class():

    seq = 'tʰ ɔ x ˈth ə r A'.split(' ')

    assert token2class(seq[0], rc('dolgo')) == 'T'
    assert token2class(seq[3], 'dolgo') == 'T'
    assert token2class(seq[-1], 'dolgo') == '0'

def test_tokens2class():

    seq = 'tʰ ɔ x ˈth ə r A ˈI'.split(' ')

    assert tokens2class(seq, 'dolgo') == list('TVKTVR00')

    assert_raises(ValueError, tokens2class, 'b  l'.split(' '), 'dolgo')

def test_prosodic_string():

    seq = 'tʰ ɔ x t ə r'.split(' ')

    assert prosodic_string(seq) == 'AXMBYN'
    
    seq = 'th o x ¹ t e'.split(' ')

    assert prosodic_string(seq) == 'AXLTBZ'

    seq = 'th o x _ th o x'.split(' ')

    assert prosodic_string(seq) == 'AXN_AXN'

    assert not prosodic_string('')

    # test for line breaks and starting with vowel
    # this is an issue in the algorithm itself!
    #seq = 'th o x t a _ th o'.split(' ')
    #assert prosodic_string(seq) == 'AXMBZ_AXLN'

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
    
    f = get_all_ngrams('ab')
    assert f == ['ab', 'a', 'b']

def test_sampa2uni():

    seq = 'tʰɔxtər'
    sampa = eval('"'+sampa2uni('t_hOxt@r')+'"')
    assert sampa == seq #or sampa2 == seq

def test_bigrams():
    
    f = bigrams('ab')
    assert f[0] == ('#','a') and f[-1] == ('b','$')

def test_trigrams():

    assert trigrams('ab')[0] == ('#', '#' ,'a') and trigrams('ab')[-1] == ('b','$', '$')

def test_fourgrams():
    
    f = fourgrams('ab')
    print(f)

    assert f[0] == ('#','#','#','a')
    assert f[-1] == ('b','$','$','$')

def test_get_n_ngrams():
    
    f = get_n_ngrams('ma',5)

    assert f[0] == ('m','a','$','$','$')

def test_pgrams():
    
    f = pgrams('ab')
    assert f[0] == ('a','X')
    assert f[-1] == ('b','N')
