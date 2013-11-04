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
    >>> tokenize(seq,model='sca')
    ['B', 'Y', 'C', 'E']

    """
    defaults = dict(
            merge_vowels = True,
            tokenizer    = False
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
            ] and not keywords['tokenizer']:
        raise ValueError('[!] Orthography {0} not available.'.format(
                    orthography)
                    )
    
    if keywords['tokenizer']:
        tfunc = lambda x: keywords['tokenizer'].graphemes(x)
    elif orthography in ['fuzzy','fuzzy_ipa']:
        tfunc = lambda x: ipa2tokens(x,**keywords)
    elif orthography in ['plain','plain_ipa']:
        tfunc = lambda x: _ipa_tokenizer.tokenize_ipa(x).split(' ')
    elif orthography in ['asjp','asjp_code']:
        tfunc = lambda x: _asjp_tokenizer.graphemes(x).split(' ')

    if not model:
        return tfunc(sequence)
    else:
        if model in rcParams:
            return tokens2class(tfunc(sequence),rcParams[model])
        else:
            return tokens2class(tfunc(sequence),model)


