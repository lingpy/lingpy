# *-* coding: utf-8 *-*
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 14:04
# modified : 2013-11-04 22:39
"""
Module provides methods and functions for dealing with linguistic sequences.
"""

__author__="Johann-Mattis List"
__date__="2013-11-04"

# import for basic import
from ..settings import rcParams
from .sound_classes import *
from .tokenizer import *
from ..util import setdefaults, data_path


_tokenizers = {}


def _get_tokenizer(name=None):
    if name not in _tokenizers:
        _tokenizers[name] = Tokenizer(data_path('orthography_profiles', name + '.prf'))
    return _tokenizers[name]


def tokenize(
        sequence,
        orthography=rcParams['basic_orthography'],
        model=False,
        **keywords):
    """
    Parameters
    ----------
    sequence : str
        The input sequence that shall be tokenized.
    orthography : {"asjp","fyzzy_ipa","plain_ipa"}(default="fuzzy_ipa")
        The name of the orthography that is assumed for the input string.
    model : {"sca","dolgo","asjp"}
        The name of the sound-class model according to whicht the string shall
        be converted.
    merge_vowels : bool (default=True)
        Specify, whether vowels shall be merged (only works with "asjp" and
        "ipa" as orthography).
    tokenizer : "str"
        Name of the orthography profile that shall be used for tokenization
        (does currently not support the conversion into sound class models).

    Examples
    --------
    
    >>> from lingpy import *
    >>> seq = "p͡fyt͡sə"
    >>> tokenize(seq)
    ['p͡f', 'y', 't͡s', 'ə']
    >>> tokenize(seq, model='sca')
    ['B', 'Y', 'C', 'E']

    """
    setdefaults(keywords, merge_vowels=True, tokenizer=False)

    if keywords['tokenizer']:
        tfunc = lambda x: keywords['tokenizer'].graphemes(x)
    elif orthography in ['fuzzy', 'fuzzy_ipa']:
        tfunc = lambda x: ipa2tokens(x, **keywords)
    elif orthography in ['plain', 'plain_ipa']:
        tfunc = lambda x: _get_tokenizer().tokenize_ipa(x).split(' ')
    elif orthography in ['asjp', 'asjp_code']:
        tfunc = lambda x: _get_tokenizer('asjp').graphemes(x).split(' ')
    else:
        raise ValueError('Orthography {0} not available.'.format(orthography))

    tokens = tfunc(sequence)
    return token2class(tokens, rcParams.get(model, model)) if model else tokens
