# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 13:27
# modified : 2013-03-04 14:28
"""
Module provides various methods for the handling of sound classes.

"""

__author__="Johann-Mattis List"
__date__="2013-03-04"


# lingpy imports
from ..data import *

def ipa2tokens(
        istring,
        diacritics = None,
        vowels = None,
        tones = None,
        combiners = '\u0361\u035c',
        breaks = '.-',
        stress = "ˈˌ'",
        merge_vowels = True
        ):
    """
    Tokenize IPA-encoded strings. 
    
    Parameters
    ----------

    seq : str
        The input sequence that shall be tokenized.
    
    diacritics : {str, None} (default=None)
        A string containing all diacritics which shall be considered in the
        respective analysis. When set to *None*, the default diacritic string
        will be used.
    
    vowels : {str, None} (default=None)
        A string containing all vowel symbols which shall be considered in the
        respective analysis. When set to *None*, the default vowel string will
        be used.

    tones : {str, None} (default=None)
        A string indicating all tone letter symbals which shall be considered
        in the respective analysis. When set to *None*, the default tone string
        will be used.

    combiners : str (default="\u0361\u035c")
        A string with characters that are used to combine two separate
        characters (compare affricates such as t͡s).

    breaks : str (default="-.")
        A string containing the characters that indicate that a new token
        starts right after them. These can be used to indicate that two
        consecutive vowels should not be treated as diphtongs or for diacritics
        that are put before the following letter.
    
    merge_vowels : bool
        Indicate, whether vowels should be merged into diphtongs
        (default=True), or whether each vowel symbol should be considered
        separately.

    Returns
    -------
    tokens : list
        A list of IPA tokens.


    Examples
    --------
    >>> from lingpy import *
    >>> myseq = 't͡sɔyɡə'
    >>> ipa2tokens(myseq)
    ['t͡s', 'ɔy', 'ɡ', 'ə']
    
    See also
    --------
    tokens2class
    class2tokens
    """

    # check for parameters
    if not vowels:
        vowels = ipa_vowels
    if not tones:
        tones = ipa_tones
    if not diacritics:
        diacritics = ipa_diacritics
    
    # create the list for the output
    out = []
    
    # set basic characteristics
    vowel = False # no vowel
    tone = False # no tone
    merge = False # no merge command
    start = True # start of unit

    for char in istring:
        # check for breaks first, since they force us to start anew
        if char in breaks:
            start = True
            vowel = False
            tone = False
            merge = False
        
        # check for combiners next
        elif char in combiners:
            out[-1] += char
            merge = True

        # check for stress
        elif char in stress:
            out += [char]
            start = True
            merge = True
            tone = False
            vowel = False
        
        # check for merge command 
        elif merge:
            out[-1] += char
            if char in vowels:
                vowel = True
            merge = False
        
        # check for diacritics
        elif char in diacritics:
            if not start:
                out[-1] += char
            else:
                out += [char]
                start = False
                merge = True
        
        # check for vowels
        elif char in vowels:
            if vowel and merge_vowels:
                out[-1] += char
            else:
                out += [char]
                vowel = True
            start = False
            tone = False
        
        # check for tones
        elif char in tones:
            if tone:
                out[-1] += char
            else:
                out += [char]
                tone = True
            start = False

        # consonants
        else:
            vowel = False
            out += [char]
            start = False
            tone = False

    return out

def tokens2class(
        tstring,
        model,
        stress = "ˈˌ'",

        ):
    """
    Convert tokenized IPA strings into their respective class strings.

    Parameters
    ----------

    tokens : list
        A list of tokens as they are returned from :py:func:`ipa2tokens`.

    model : :py:class:`~lingpy.data.model.Model`
        A :py:class:`~lingpy.data.model.Model` object.

    Returns
    -------

    classes : string
        A sound-class representation of the tokenized IPA string.

    Examples
    --------
    >>> from lingpy import *
    >>> tokens = ipa2tokens('t͡sɔyɡə')
    >>> classes = tokens2class(tokens,sca)
    >>> print(classes)
    CUKE

    See also
    --------
    ipa2tokens
    class2tokens

    """
    out = []
    for token in tstring:
        try:
            out.append(model.converter[token])
        except:
            try:
                out.append(model.converter[token[0]])
            except:
                # check for stressed syllables
                if token[0] in stress:
                    try:
                        out.append(model.converter[token[1:]])
                    except:
                        out.append(model.converter[token[1]])
                else:
                    # new character for missing data and spurious items
                    out.append('0')

    return out

def prosodic_string(
        tstring
        ):
    """
    Create a prosodic string of the sonority profile of a sequence.

    Paremeters
    ----------

    seq : list
        A list of integers indicating the sonority of the tokens of the
        underlying sequence.

    Returns
    -------
    prostring : string
        A prosodic string corresponding to the sonority profile of the
        underlying sequence.

    See also:
    ---------
    
    prosodic weights

    Notes
    -----
    
    A prosodic string is a sequence of specific characters which indicating
    their resprective prosodic context (see :evobib:`List2012` or
    :evobib:`List2012a` for a detailed description).


    Examples
    --------
    >>> profile = [int(i) for i in tokens2class(ipa2tokens('t͡sɔyɡə'),art)]
    >>> prosodic_string(profile)
    '#vC>'

    """
    # get the sonority profile
    sstring = [int(t) for t in tokens2class(tstring,art)]

    # create a trigram-representation of the string
    trigrams = zip(
            [0] + sstring,
            sstring,
            sstring[1:] + [0]
            )
    # in this updated version, we will produce more meta-data, in order to
    # avoid being restricted to one single model
    
    # set up some important values
    syl = 1      # syllable counter
    pstring = [] # prosodic string
    
    for a,b,c in trigrams:
        
        # check for vowel, vowels are always peaks
        if b == 7:
            pstring += ['^']

        # check for peak (consonants)
        elif a < b > c:
            pstring += ['^']

        # check for descending position
        elif a >= b >= c:
            pstring += ['>']

        # ascending #1
        elif b < c or a > b <= c or a < b <= c:
            pstring += ['<']

        # dummy for other stuff
        else:
            pstring += ['?']

    
    return pstring
            
