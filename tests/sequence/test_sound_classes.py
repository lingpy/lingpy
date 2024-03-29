"""
Tests for the SCA module
"""
import pytest

from lingpy import rc, csv2list
from lingpy.sequence.sound_classes import ipa2tokens, token2class, \
    tokens2class, prosodic_string, prosodic_weights, class2tokens, pid, \
    check_tokens, sampa2uni, pgrams, syllabify, tokens2morphemes, ono_parse, \
    clean_string, _get_brackets, codepoint


def test_ipa2tokens(test_data):
    seq = 'ˈtʲʰoɔːix_tərp͡f¹¹'
    assert len(ipa2tokens(seq)) != len(list(seq))

    seq = 'ʰto͡i'
    assert len(ipa2tokens(seq)) == 2

    with pytest.raises(ValueError):
        seq = ['t͡s', 'ɔ', 'y', 'ɡ', 'ə']
        ipa2tokens(seq)

    with pytest.raises(ValueError):
        seq = 'th o x t a'
        ipa2tokens(seq)

    # now check with all possible data we have so far, but only on cases
    # where tokenization doesn't require the merge_vowels = False flag
    tokens = csv2list(str(test_data / 'test_tokenization.tsv'))

    for a, b in tokens:
        tks = ' '.join(ipa2tokens(a))
        # we check for two variants, since we don't know whether vowels are
        # merged or not in the test data
        assert tks == b

    # now test on smaller set with unmerged vowels
    tokens = csv2list(str(test_data / 'test_tokenization_mv.tsv'))

    for a, b in tokens:
        tks = ' '.join(ipa2tokens(a, merge_vowels=False,
                                  merge_geminates=False))
        # we check for two variants, since we don't know whether vowels are
        # merged or not in the test data
        assert tks == b

    tokens = csv2list(str(test_data / 'test_tokenization_nasals.tsv'))
    for a, b in tokens:
        tks = ' '.join(ipa2tokens(a, merge_vowels=True,
                                  merge_geminates=True, expand_nasals=True,
                                  semi_diacritics='h'))
        assert tks == b


def test_token2class():
    seq = 'tʰ ɔ x ˈth ə r A'.split(' ')

    assert token2class(seq[0], rc('dolgo')) == 'T'
    assert token2class(seq[3], 'dolgo') == 'T'
    assert token2class(seq[-1], 'dolgo') == '0'
    assert token2class('', 'dolgo') == '0'


def test_tokens2class():
    seq = 'tʰ ɔ x ˈth ə r A ˈI ʲ'.split(' ')
    seq2 = 'th o ?/x a'.split(' ')
    seq3 = 'th o ?/ a'.split(' ')

    assert tokens2class(seq, 'dolgo') == list('TVKTVR000')
    assert tokens2class(seq2, 'cv', cldf=False)[2] == '0'
    assert tokens2class(seq2, 'cv')[2] == 'C'
    assert tokens2class(seq3, 'cv', cldf=True)[2] == '0'

    with pytest.raises(ValueError):
        tokens2class(['A'], 'dolgo')
    with pytest.raises(ValueError):
        tokens2class('bla', 'sca')


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
    # seq = 'th o x t a _ th o'.split(' ')
    # assert prosodic_string(seq) == 'AXMBZ_AXLN'


def test_prosodic_weights():
    seq = 'tʰ ɔ x t ə r'.split(' ')

    assert prosodic_weights(prosodic_string(seq))[0] == 2
    assert prosodic_weights(prosodic_string(seq))[-1] == 0.8


def test_class2tokens():
    classes = 'T-VKTV-R'
    tokens = 'tʰ ɔ x t ə r'.split(' ')

    out = class2tokens(classes, tokens)
    _ = class2tokens([['T'], ['-VKTV-'], ['R']], 'th o x t e r'.split(),
                     local=True)

    assert out[1] == '-' and out[-2] == '-'


def test_pid():
    assert pid('mattis', 'maTTIs', 1) == 0.5
    assert pid('mattis', 'maTTIs', 2) == 0.5
    assert pid('mattis', 'maTTIs', 3) == 0.5
    assert pid('mattis', 'maTTIs', 4) == 0.5
    assert pid('m-', '-m', mode=1) == 0
    assert pid('m-', '-m', mode=2) == 0
    assert pid('m', '-', mode=3) == 0
    assert pid('m-', '-m', mode=4) == 0


def test_check_tokens():
    assert check_tokens('th o x T e r'.split(' '))[0] == (3, 'T')


def test_sampa2uni():
    seq = 'tʰɔxtər'
    sampa = eval('"' + sampa2uni('t_hOxt@r') + '"')
    assert sampa == seq  # or sampa2 == seq


def test_pgrams():
    f = pgrams('ab')
    assert f[0] == ('a', 'X')
    assert f[-1] == ('b', 'N')


def test_syllabify():
    seq1 = "t i a o ¹ b u ² d a o"
    seq2 = "jabloko"
    seq3 = "jabəlko"
    seq4 = "j a b u - k o"
    seq5 = "ma⁵io"

    with pytest.raises(ValueError):
        syllabify(seq1, output="test")

    assert syllabify(seq1, output="flat").count(rc('morpheme_separator')) == 2
    assert syllabify(seq2, output="breakpoints")[0] == (0, 2)
    assert syllabify(seq3, output="nested")[1] == list("bəl")
    assert syllabify(seq4, output="nested")[1] == list("bu-")
    assert ''.join(syllabify(seq5, sep="+")).split('+')[-1] == 'io'


def test_tokens2morphemes():
    seq1 = "t i a o ¹ b u ² d a o".split(' ')
    seq2 = "t i a o ¹ + b u ² # d a o".split(' ')
    seq3 = "t i a o ¹ b u _ d a o".split(' ')
    seq4 = "t i a o murks w a o".split(' ')

    assert len(tokens2morphemes(seq1)) == 3
    assert len(tokens2morphemes(seq2)) == 3
    assert len(tokens2morphemes(seq3)) == 2
    assert len(tokens2morphemes(seq4, sep='murks')) == 2
    assert len(tokens2morphemes(seq1, split_on_tones=False)) == 1
    with pytest.raises(ValueError):
        tokens2morphemes("t i a o")
    assert len(tokens2morphemes(list("b++t"))) == 2


def test_onoparse():
    seq = "a k e r ts a n".split(' ')
    out1 = ono_parse(seq, output='pprint')
    out2 = ono_parse(seq, output='prostring')
    out3 = ono_parse(seq)

    assert isinstance(out1, str)
    assert out3[0] == ('-', '#')
    assert out2 == 'VCvcC>$'


def test_clean_string():
    seq1 = 'this (is an error)'
    seq2 = 'feature/vector'
    seq3 = 'ta:tata'
    seq4 = 'what (the) hack [this is]'
    seq5 = 't a t'

    _get_brackets('A')

    assert clean_string(seq1)[0] == 'th i s'
    assert clean_string(seq2)[1] == 'v e c t o r'
    assert clean_string(seq3)[0] == 't a: t a t a'
    assert clean_string(seq4)[0] == 'wh a t _ h a c k'
    assert clean_string(seq5, segmentized=True)[0] == 't a t'
    assert clean_string('a(a', ignore_brackets=False)[0] == 'a ( a'
    assert clean_string('a/a', split_entries=False)[0] == 'a / a'
    assert clean_string('aa', preparse=[('a', 'b')])[0] == 'bb'
    assert clean_string('bb', merge_geminates=False)[0] == 'b b'
    assert clean_string('bb', rules={"b": "cd"},
                        merge_geminates=False)[0] == "cd cd"

def test_codepoint():
    assert codepoint('á') == 'U+00e1'
