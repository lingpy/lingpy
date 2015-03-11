# *-* coding: utf-8 *-*
# These lines are for py2-py3 compatibility
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

"""
Module provides various methods for the handling of sound classes.
"""

# lingpy imports
import re
from ..settings import rcParams
from ..data.ipa.sampa import reXS,xs

def ipa2tokens(
        istring,
        **keywords
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
    
    merge_vowels : bool (default=False)
        Indicate, whether vowels should be merged into diphtongs
        (default=True), or whether each vowel symbol should be considered
        separately.

    merge_geminates : bool (default=False)
        Indicate, whether identical symbols should be merged into one token, or
        rather be kept separate.

    excand_nasals : bool (default=False)

    semi_diacritics: str (default='')
        Indicate which symbols shall be treated as "semi-diacritics", that is,
        as symbols which can occur on their own, but which eventually, when
        preceded by a consonant, will form clusters with it. If you want to
        disable this features, just set the keyword to an empty string.
    clean_string : bool (default=False)
        Conduct a rough string-cleaning strategy by which all items between
        brackets are removed along with the brackets, and 

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
    # go for defaults
    kw = dict(
            breaks = rcParams['breaks'],
            combiners = rcParams['combiners'],
            diacritics = rcParams['diacritics'],
            expand_nasals = False,
            merge_geminates = True,
            merge_vowels = rcParams['merge_vowels'],
            semi_diacritics = '',
            stress = rcParams['stress'],
            tones = rcParams['tones'],
            vowels = rcParams['vowels'],
            clean_sequence = False # add this later, not today XXX
            )
    kw.update(keywords)
    
    if kw['clean_sequence']:
        raise ValueError("This part has not yet been implemented!")

    # check for pre-tokenized strings
    if ' ' in istring:
        out = istring.split(' ')
        if istring.startswith('#'):
            return out[1:-1]
        else:
            return out
    
    # create the list for the output
    out = []
    
    # set basic characteristics
    vowel = False # no vowel
    tone = False # no tone
    merge = False # no merge command
    start = True # start of unit
    nasal = False
    
    # define nasals and nasal chars and semi_diacritics
    nasals = 'ãũẽĩõ'
    nasal_char = "\u0303"
    semi_diacritics = kw['semi_diacritics']
    nogos = "_◦"


    for char in istring:
        # check for nasal stack and vowel environment
        if nasal:
            if char not in kw['vowels'] and char not in kw['diacritics']:
                out += [rcParams['nasal_placeholder']]
                nasal = False

        # check for breaks first, since they force us to start anew
        if char in kw['breaks']:
            start = True
            vowel = False
            tone = False
            merge = False
        
        # check for combiners next
        elif char in kw['combiners']:
            out[-1] += char
            merge = True

        # check for stress
        elif char in kw['stress']:
            out += [char]
            # XXX be careful about the removement of the start-flag here, but it
            # XXX seems to make sense so far!
            merge = True
            tone = False
            vowel = False
            start = False
        
        # check for merge command 
        elif merge:
            out[-1] += char
            if char in kw['vowels']:
                vowel = True
            merge = False

        # check for nasals in NFC normalization and non-normalizable nasals
        elif kw['expand_nasals'] and char == nasal_char and vowel:
            out[-1] += char
            start = False
            nasal = True

        # check for weak diacritics
        elif char in semi_diacritics and not start and not vowel and not tone and out[-1] not in nogos:
            out[-1] += char
        
        # check for diacritics
        elif char in kw['diacritics']:
            if not start:
                out[-1] += char
            else:
                out += [char]
                start = False
                merge = True
        
        # check for vowels
        elif char in kw['vowels']:
            if vowel and kw['merge_vowels']:
                out[-1] += char
            else:
                out += [char]
                vowel = True
            start = False
            tone = False

            if kw['expand_nasals'] and char in nasals:
                nasal = True
        
        # check for tones
        elif char in kw['tones']:
            vowel = False
            if tone:
                out[-1] += char
            else:
                out += [char]
                tone = True
            start = False

        # consonants
        else:
            vowel = False
            tone = False
            out += [char]
            start = False
            tone = False

    if nasal:
        out += [rcParams['nasal_placeholder']]

    if kw['merge_geminates']:
        new_out = [out[0]]
        for i in range(len(out) -1):
            outA = out[i]
            outB = out[i+1]
            if outA == outB:
                new_out[-1] += outB
            else:
                new_out += [outB]
        return new_out
    
    return out

def syllabify(seq, alignment=False, breakpoints=False):
    """
    Carry out a simple syllabification of a sequence, using sonority as a proxy.

    Note
    ----
    When analyzing the sequence, we start a new syllable in all cases where we
    reach a deepest point.

    Returns
    -------
    syllable : list
        A nested list, representing the structure of the syllable.
    """
    # check for alignment keyword
    if alignment:

        if ' ' in seq:
            seq = seq.split(' ')
        alm = [x for x in seq]
        
        seq = [s for s in seq if s != '-']
    
    # we assume we are dealing with tokens if the syllable is a list.
    if type(seq) == list:
        listed_seq = [s for s in seq]
    elif ' ' in seq:
        listed_seq = seq.split(' ')
    else:
        listed_seq = ipa2tokens(seq)

    profile = [0] + [int(i) for i in tokens2class(listed_seq, 'art')] + [0]
    
    new_syl = False
    breaks = []
    for i in range(1,len(profile)-1):
        # get the pro-tokens
        p1,p2,p3 = profile[i-1],profile[i],profile[i+1]

        # get the char
        char = listed_seq[i-1]
        
        # simple rule: we start a new syllable, if p2 is smaller or equal to p1 and p3
        # is larger than p2
        if p1 >= p2 < p3:
            new_syl = True

        # get the char before, after
        if new_syl:
            breaks += [rcParams['morpheme_separator'],char]
            new_syl = False
        else:
            breaks += [char]
    
    # if the alignment option is chosen, we need to re-assign the gaps to the
    # current form
    if alignment:
        out = []
        idxA,idxB = 0,0
        bpoints = []
        lastbp = 0
        
        while idxB < len(alm) and idxA < len(breaks):
            charA = breaks[idxA]
            charB = alm[idxB]

            if charA == charB:
                idxA += 1
                idxB += 1
                out += [charA]

            elif charB == '-':
                idxB += 1
                out += [charB]

            elif charA == rcParams['morpheme_separator']:
                idxA += 1
                out += [charA]
                bpoints += [(lastbp,len(out)-1)]
                lastbp = len(out)-1

        bpoints += [(lastbp,len(out))]

        if not breakpoints:
            return out
        else:
            return bpoints

        


    return breaks

def asjp2tokens(
        seq,
        merge_vowels = True
        ):
    """
    """
    
    tokens = ' '.join(
                ipa2tokens(
                    seq,
                    diacritics = '*$~"',
                    vowels = 'aeiouE3',
                    tones = '',
                    combiners = '',
                    merge_vowels = merge_vowels
                    )
                )
    tokens = re.sub(r'([^ ]) ([^ ])~',r'\1\2~',tokens)
    tokens = re.sub(r'([^ ]) ([^ ]) ([^ ])\$',r'\1\2\3$',tokens)

    return tokens.split(' ')

def token2class(
        token,
        model,
        **keywords
        ):
    """
    Convert a single token into a sound-class.

    tokens : str
        A token (IPA-string).

    model : :py:class:`~lingpy.data.model.Model`
        A :py:class:`~lingpy.data.model.Model` object.

    Returns
    -------
    c : str
        The corresponding sound-class value.

    """
    kw = dict(
            stress = rcParams['stress']
            )
    kw.update(keywords)

    # check whether model is passed as real model or as string
    if str(model) == model:
        model = rcParams[model]

    try:
        return model[token]
    except:
        try:
            return model[token[0]]
        except:
            try:
                # check for stress
                return model[token[1:]]
            except:
                try:
                    return model[token[1]]
                except:
                    return '0'

def tokens2class(
        tstring,
        model,
        **keywords
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
    >>> classes = tokens2class(tokens,'sca')
    >>> print(classes)
    CUKE

    See also
    --------
    ipa2tokens
    class2tokens

    """
    kw = dict(
            stress = rcParams['stress']
            )
    kw.update(keywords)
    
    # check whether model is passed as real model or as string
    if str(model) == model:
        model = rcParams[model]
    
    out = []
    for token in tstring:
        try:
            out.append(model[token])
        except:
            try:
                out.append(model[token[0]])
            except:
                # check for stressed syllables
                if len(token) > 0:
                    if token[0] in kw['stress']:
                        try:
                            out.append(model[token[1:]])
                        except:
                            try:
                                out.append(model[token[1]])
                            except:
                                # new character for missing data and spurious items
                                out.append('0')
                    else:
                        # new character for missing data and spurious items
                        out.append('0')
                else:
                    raise ValueError("[!] string '{0}' is erroneously coded!".format(tstring))

    return out

def prosodic_string(
        string,
        _output = True,
        **keywords
        ):
    """
    Create a prosodic string of the sonority profile of a sequence.

    Parameters
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
    In contrast to the previous model, the current implementation allows for a
    more fine-graded distinction between different prosodic segments. The
    current scheme distinguishes 9 prosodic positions:
    
    * ``A``: sequence-initial consonant
    * ``B``: syllable-initial, non-sequence initial consonant in a context of
      ascending sonority
    * ``C``: non-syllable, non-initial consonant in ascending sonority context
    * ``L``: non-syllable-final consonant in descending environment
    * ``M``: syllable-final consonant in descending environment
    * ``N``: word-final consonant
    * ``X``: first vowel in a word
    * ``Y``: non-final vowel in a word
    * ``Z``: vowel occuring in the last position of a word

    Examples
    --------
    >>> prosodic_string(ipa2tokens('t͡sɔyɡə')
    'AXBZ'

    """
    defaults = dict(
            stress = rcParams['stress']
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # check for empty string passed
    if not string:
        return ''

    # check for the right string
    if not isinstance(string[0], int):
        
        # get the sonority profile
        sstring = [9] + [int(t) for t in tokens2class(
            string,
            rcParams['art'],
            **keywords)
            ] + [9]
    else:
        sstring = [9] + string + [9]

    # check for multiple strings in string
    if 9 in sstring[1:-1]:

        # break the string into pieces
        nstrings = [[]]
        switch = False
        for i in sstring[1:-1]:
            if i != 9:
                nstrings[-1] += [i]
            else:
                nstrings += [[]]
        
        # return the prostrings of the pieces recursively, note that the
        # additional check whether x is True is necessitated by the fact that
        # often errors occur in the coding, i.e. strings are given 
        return '_'.join([prosodic_string(x,_output) for x in nstrings])
    
    # create the output values
    pstring = ''
    first = True # stores whether first syllable is currently being processed

    # start iterating over relevant parts of the string
    for i in range(1,len(sstring)-1):

        # get all three values
        a,b,c = sstring[i-1],sstring[i],sstring[i+1]

        # check for vowel
        if b == 7:
            
            # check for first
            if first:
                pstring += 'X'
                first = False

            # check for last
            elif c == 9:
                pstring += 'Z'

            else:
                pstring += 'Y'
        
        # check for tones
        elif b == 8:
            pstring += 'T'

        # check for descending position
        elif a >= b >= c or c == 8:

            # check for word final consonant
            if c == 9 and b != 7:
                pstring += 'N'
            
            # check for word final vowel
            elif c == 9 and b == 7:
                pstring += 'Z'
            
            else:
                if first:
                    first = False
                    pstring += 'A'
                else:
                    pstring += 'L'

        # ascending
        elif b < c or a > b <= c or a < b <= c:

            # check for syllable first
            if a == 9:
                pstring += 'A'
            elif a >= b:
                if c == 9:
                    pstring += 'N'
                else:
                    if pstring[-1] != 'A':
                        pstring = pstring[:-1]+pstring[-1].replace('L','M')+'B'
                    else:
                        pstring += 'C'
            else:
                pstring += 'C'
        
        # consonant peak
        elif a < b > c:
            if first:
                pstring += 'X'
                first = False
            else:
                pstring += 'Y'

        # dummy for other stuff
        else:
            print("[i] Warning, condition not met in conversion.")
            print(sstring)
            print(pstring)
            print(a,b,c)
            input("[i] Press the Any-Key to carry on.")
            pstring += '?'
    
    if _output == 'cv':
        conv = {
                "A":"C",
                "B":"C",
                "C":"C",
                "M":"C",
                "L":"C",
                "N":"C",
                "X":"V",
                "Y":"V",
                "Z":"V",
                "T":"T",
                "_":"_",
                }
        return ''.join([conv[x] for x in pstring])

    elif _output == 'CcV':
        conv = {
                "A":"C",
                "B":"C",
                "C":"C",
                "M":"c",
                "L":"c",
                "N":"c",
                "X":"V",
                "Y":"V",
                "Z":"V",
                "T":"T",
                "_":"_",
                }
        return ''.join([conv[x] for x in pstring])

    elif _output:
        return pstring

    else:
        conv = {
                "A":"#",
                "B":"C",
                "C":"C",
                "M":"c",
                "L":"c",
                "N":"$",
                "X":"V",
                "Y":"v",
                "Z":">"
                }
        return ''.join([conv[x] for x in pstring])

def prosodic_weights(
        prostring,
        _transform = {}
        ):
    """
    Calculate prosodic weights for each position of a sequence.
    
    Parameters
    ----------

    prostring : string
        A prosodic string as it is returned by :py:func:`prosodic_string`.
    _transform : dict
        A dictionary that determines how prosodic strings should be transformed
        into prosodic weights. Use this dictionary to adjust the prosodic
        strings to your own user-defined prosodic weight schema.

    Returns
    -------
    weights : list
        A list of floats reflecting the modification of the weight for each position.

    Notes
    -----

    Prosodic weights are specific scaling factors which decrease or increase
    the gap score of a given segment in alignment analyses (see :evobib:`List2012` or
    :evobib:`List2012a` for a detailed description).

    Examples
    --------
    >>> from lingpy import *
    >>> prostring = '#vC>'
    >>> prosodic_weights(prostring)
    [2.0, 1.3, 1.5, 0.7]
    
    See also
    --------
    prosodic_string

    """
    # check for transformer
    if _transform:
        transform = _transform
    
    # default scale for tonal languages
    elif 'T' in prostring:
        transform = {
                '#' : 1.6,
                'V' : 3.0,
                'C' : 1.2,
                'c' : 1.1,
                'v' : 3.0, # raise the cost for the gapping of vowels
                '<' : 0.8, 
                '$' : 0.5, 
                '>' : 0.7, 
                'T' : 1.0,
                '_' : 0.0,

                # new values for alternative prostrings
                'A' : 1.6,  # initial
                'B' : 1.3, # syllable-initial
                'C' : 1.2,  # ascending

                'L' : 1.1,  # descending
                'M' : 1.1,  # syllable-descending
                'N' : 0.5,  # final
                
                'X' : 3.0,  # vowel in initial syllable
                'Y' : 3.0,  # vowel in non-final syllable
                'Z' : 0.7,  # vowel in final syllable
                'T' : 1.0,  # Tone
                '_' : 0.0   # break character
                }
    # default scale for other languages
    else:
        transform = {
                '#' : 2.0,
                'V' : 1.5,
                'C' : 1.5,
                'c' : 1.1,
                'v' : 1.3,
                '<' : 0.8, 
                '$' : 0.8, 
                '>' : 0.7, 
                'T' : 0.0,
                '_' : 0.0,

                # new values for alternative prostrings
                'A' : 2.0,  # initial
                'B' : 1.75, # syllable-initial
                'C' : 1.5,  # ascending

                'L' : 1.1,  # descending
                'M' : 1.1,  # syllable-descending
                'N' : 0.8,  # final
                
                'X' : 1.5,  # vowel in initial syllable
                'Y' : 1.3,  # vowel in non-final syllable
                'Z' : 0.8,  # vowel in final syllable
                'T' : 0.0,  # Tone
                '_' : 0.0   # break character

                }
    
    out = [transform[i] for i in prostring]

    return out

def class2tokens(
        tokens,
        classes,
        gap_char = '-',
        local = False
        ):
    """
    Turn aligned sound-class sequences into an aligned sequences of IPA tokens.

    Parameters
    ----------

    tokens : list
        The list of tokens corresponding to the unaligned IPA string.

    classes : string or list
        The aligned class string.

    gap_char : string (default="-")
        The character which indicates gaps in the output string.
    
    local : bool (default=False)
        If set to *True* a local alignment with prefix and suffix can be
        converted.
    
    Returns
    -------
    alignment : list
        A list of tokens with gaps at the positions where they occured in the
        alignment of the class string.

    See also
    --------
    ipa2tokens
    tokens2class
    
    Examples
    --------
    >>> from lingpy import *
    >>> tokens = ipa2tokens('t͡sɔyɡə')
    >>> aligned_sequence = 'CU-KE'
    >>> print ', '.join(class2tokens(tokens,aligned_sequence))
    t͡s, ɔy, -, ɡ, ə

    """
    if not local:
        out = [t for t in tokens]
        for i in range(len(classes)):
            if classes[i] in '-X':
                out.insert(i,gap_char)
    else:        
        # get the length of the prefix
        prefix = len(classes[0])

        # get the suffix 
        suffix = -len(classes[2])

        if suffix == 0:
            suffix = None

        # get the substring
        out = [t for t in tokens[prefix:suffix]]

        # start the loop
        for i in range(len(classes[1])):
            if classes[1][i] in '-X':
                out.insert(i,gap_char)            
    return out

def pid(
        almA,
        almB,
        mode=2
        ):
    """
    Calculate the Percentage Identity (PID) score for aligned sequence pairs.

    Parameters
    ----------

    almA, almB : string or list
        The aligned sequences which can be either a string or a list.

    mode : { 1, 2, 3, 4, 5 }
        Indicate which of the four possible PID scores described in :evobib:`Raghava2006`
        should be calculated, the fifth possibility is added for linguistic
        purposes:
        
        1. identical positions / (aligned positions + internal gap positions),
        
        2. identical positions / aligned positions,
        
        3. identical positions / shortest sequence, or
        
        4. identical positions / shortest sequence (including internal gap
           pos.)

        5. identical positions / (aligned positions + 2 * number of gaps)  

    Returns
    -------

    score : float
        The PID score of the given alignment as a floating point number between
        0 and 1.

    Notes
    -----
    
    The PID score is a common measure for the diversity of a given alignment.
    The implementation employed by LingPy follows the description of
    :evobib:`Raghava2006` where four different variants of PID scores are
    distinguished. Essentially, the PID score is based on the comparison of
    identical residue pairs with the total number of residue pairs in a given
    alignment.  

    Examples
    --------
    Load an alignment from the test suite.

    >>> from lingpy import *
    >>> pairs = PSA(get_file('test.psa'))

    Extract the alignments of the first aligned sequence pair.

    >>> almA,almB,score = pairs.alignments[0]

    Calculate the PID score of the alignment.

    >>> pid(almA,almB)
    0.44444444444444442

    See also
    --------
    lingpy.compare.Multiple.get_pid

    .. todo:: change debug for ZeroDivisionError

    """

    zipped = zip(almA,almB)
    idn_pos = 0
    int_gps = 0
    aln_pos = 0

    for charA,charB in zipped:
        tmp = [charA,charB].count('-')
        if tmp == 1:
            int_gps += 1
        elif tmp == 0 and charA == charB:
            idn_pos += 1
            aln_pos += 1
        elif tmp == 0:
            aln_pos += 1

    if mode == 2:
        try:
            return idn_pos / (aln_pos + int_gps)
        except ZeroDivisionError:
            #print('\t'.join(almA))
            #print('\t'.join(almB))
            #print('-----')
            return 0

    elif mode == 1: 
        try:
            return idn_pos / aln_pos
        except ZeroDivisionError:
            #print('\t'.join(almA))
            #print('\t'.join(almB))
            #print('-----')
            return 0

    elif mode == 3:
        srt_seq = min(
                len(
                    [i for i in almA if i != '-']
                    ),
                len(
                    [i for i in almB if i != '-']
                    )
                )
        try:
            return idn_pos / srt_seq
        except ZeroDivisionError:
            #print('\t'.join(strA))
            #print('\t'.join(strB))
            #print('-----')
            return 0

    elif mode == 4:
        srt_seq = min(
                len(
                    ''.join([i[0] for i in almA]).strip('-')
                    ),
                len(
                    ''.join([i[0] for i in almB]).strip('-')
                    )
                )
        try:
            return idn_pos / srt_seq
        except ZeroDivisionError:
            print('\t'.join(almA))
            print('\t'.join(almB))
            print('-----')
            return 0

    elif mode == 5:
        
        return idn_pos / len(almA)

def check_tokens(tokens,**keywords):
    """
    Function checks whether tokens are given in a consistent input format.
    """
    defaults = dict(
            stress = rcParams['stress']
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    errors = []

    for i,token in enumerate(tokens):
        
        # check for conversion within the rcParams['art']iculation-model
        try:
            rcParams['art'].converter[token]
        except KeyError:
            try:
                rcParams['art'].converter[token[0]]
            except KeyError:
                if token[0] in keywords['stress']:
                    try:
                        rcParams['art'].converter[token[1]]
                    except KeyError:
                        errors.append((i,token))
                else:
                    errors.append((i,token))
    
    return errors

def get_all_ngrams(sequence):
    """
    Function returns all possible n-grams of a given sequence.

    Parameters
    ----------
    sequence : list or str
        The sequence that shall be converted into it's ngram-representation.

    Returns
    -------
    out : list
        A list of all ngrams of the input word, sorted in decreasing order of
        length.

    Examples
    --------
    >>> get_all_ngrams('abcde')
    ['abcde', 'bcde', 'abcd', 'cde', 'abc', 'bcd', 'ab', 'de', 'cd', 'bc', 'a', 'e', 'b', 'd', 'c']

    """
    
    # get the length of the word
    l = len(sequence)

    # determine the starting point 
    i = 0

    # define the output list
    out = []

    # start the while loop
    while i != l and i < l:
        # copy the sequence
        new_sequence = sequence[i:l]

        # append the sequence to the output list
        out += [new_sequence]

        # loop over the new sequence
        for j in range(1,len(new_sequence)):
            out += [new_sequence[:j]]
            out += [new_sequence[j:]]

        # increment i and decrement l
        i += 1
        l -= 1

    return sorted(out,key=lambda x:len(x),reverse=True)

def sampa2uni(seq):
    """
    Convert sequence in IPA-sampa-format to IPA-unicode.

    Notes
    -----
    This function is based on code taken from Peter Kleiweg
    (http://www.let.rug.nl/~kleiweg/L04/devel/python/xsampa.html). 

    """

    result = ''
    tokens = reXS.findall(seq)
    for tok, err in tokens:
        assert not err and tokens
        result += xs[tok]

    return result

def bigrams(sequence):
    """
    Convert a given sequence into a sequence of bigrams.
    """
    if ' ' in sequence:
        seq = sequence.split(' ')
    else:
        seq = list(sequence)

    return list(zip(['#']+seq,seq+['$']))

def trigrams(sequence):
    """
    Convert a given sequence into a sequence of trigrams.
    """
    if ' ' in sequence:
        seq = sequence.split(' ')
    else:
        seq = list(sequence)

    return list(
            zip(
                ['#','#']+seq,
                ['#']+seq+['$'],
                seq+['$','$']
                )
            )

def fourgrams(sequence):
    """
    Convert a given sequence into a sequence of trigrams.
    """
    if ' ' in sequence:
        seq = sequence.split(' ')
    else:
        seq = list(sequence)

    return list(
            zip(
                ['#','#','#']+seq,
                ['#','#']+seq+['$'],
                ['#']+seq+['$','$'],
                seq+['$','$','$']
                )
            )

def get_n_ngrams(sequence, ngram=4):
    """
    convert a given sequence into a sequence of ngrams.
    """

    if ' ' in sequence:
        seq = sequence.split(' ')
    else:
        seq = list(sequence)

    
    tobezipped = []
    for i in range(ngram):
        prefix = (ngram-i-1) * ['#']
        postfix = i * ['$']
        tobezipped += [tuple(prefix + seq + postfix)]
    
    return list(zip(*tobezipped))[ngram-1:]

def pgrams(sequence, **keywords):
    """
    Convert a given sequence into bigrams consisting of prosodic string symbols and the tokens of the original sequence.
    """
    if ' ' in sequence:
        seq = sequence.split(' ')
    elif isinstance(sequence, str):
        seq = ipa2tokens(sequence)
    else:
        seq = list(sequence)

    return list(
            zip(
                seq,
                prosodic_string(seq, **keywords)
                )
            )
