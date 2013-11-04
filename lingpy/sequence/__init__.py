# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 14:04
# modified : 2013-03-04 14:04
"""
Module provides methods and functions for dealing with linguistic sequences.
"""

__author__="Johann-Mattis List"
__date__="2013-03-04"

# import for basic import
from ..settings import rcParams
from .sound_classes import *
from .tokenizer import *

# define a static tokenizer for the tokenize-function
_ipa_tokenizer = Tokenizer()
_asjp_tokenizer = Tokenizer('asjp')

# define function for easy tokenization
def tokenize(
        sequence,
        orthography=rcParams['basic_orthography'],
        model=False,
        **keywords
        ):
    defaults = dict(
            merge_vowels = True
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # check for appropriate input
    if orthography not in [
            'fuzzy_ipa',
            'plain_ipa',
            'fuzzy',
            'plain',
            'asjp'
            ]:
        raise ValueError('[!] Orthography {0} not available.'.format(
                    orthography)
                    )

    # return split of already tokenized stuff if string is already tokenized
    if ' ' in sequence:
        tokens = sequence.split(' ')
    elif orthography in ['fuzzy','fuzzy_ipa']:
        tokens = ipa2tokens(sequence)
    elif orthography in ['plain','plain_ipa']:
        tokens = _ipa_tokenizer.graphemes(sequence)
    elif orthography in ['asjp']:
        tokens = _asjp_tokenizer(sequence)
    
    # check for existing model
    if not model:
        return tokens
    elif model in rcParams:
        return tokens2class(tokens,rcParams[model])
    else:
        return tokens2class(tokens,model)


