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
from ..data.ipa.sampa import reXS,xs

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
        string,
        output = 'tuples'
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
    # check for the right string
    if type(string[0]) != int:
        
        # get the sonority profile
        sstring = [9] + [int(t) for t in tokens2class(string,art)] + [9]
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
        if output in ['prostring','p','pstring']:
            return '_'.join([prosodic_string(x,'p') for x in nstrings])
        else:
            idx = 0
            outstrings,outsyllables = [],[]
            for n in nstrings:
                pstring,syllables = prosodic_string(n,'t')
                outstrings += [pstring]
                outsyllables += [s+idx for s in syllables]
                idx += len(pstring)
            return '_'.join(outstrings),outsyllables
    
    # create the output values
    pstring = ''
    s = True
    syllables = []

    # start iterating over relevant parts of the string

    for i in range(1,len(sstring)-1):

        # get all three values
        a,b,c = sstring[i-1],sstring[i],sstring[i+1]

        # check for start
        if a == 9:
            s = True

        # check for vowel
        if b == 7:
            pstring += 'v'
        
        # check for tones
        elif b == 8:
            pstring += 'T'

        # check for peak (consonants)
        elif a < b > c:
            pstring += 'C'

        # check for descending position
        elif a >= b >= c or c == 9:
            pstring += 'c'

        # ascending #1
        elif b < c or a > b <= c or a < b <= c:
            pstring += 'C'

            if a >= b:
                s = True

        # dummy for other stuff
        else:
            pstring += '?'

        if s:
            syllables += [i-1]
            s = False

    if output in ["tuples","tuple",'t']:
        return pstring,syllables
    
    if output in ['prostring','pstring','p']:
        cA = {
                'c':'#',
                'C':'#',
                'v':'V'
                }
        cB = {
                'c':'$',
                'C':'$',
                'v':'>'
                }

        if len(syllables) == 1:
            out = cA[pstring[0]] + pstring[1:-1] + cB[pstring[-1]]
            return out.replace('v','V')
        elif len(syllables) == 2:
            return cA[pstring[0]] + pstring[1:syllables[1]].replace('v','V') +\
                    pstring[syllables[1]:-1] + cB[pstring[-1]]
        else:
            return cA[pstring[0]] + pstring[1:syllables[1]].replace('v','V') +\
                    pstring[syllables[1]:syllables[-1]] + pstring[syllables[-1]:-1] +\
                    cB[pstring[-1]]

def prosodic_weights(
        prostring,
        scale = None,
        factor = 0.3
        ):
    """
    Calculate prosodic weights for each position of a sequence.
    
    Parameters
    ----------

    prostring : string
        A prosodic string as it is returned by :py:func:`prosodic_string`.

    scale : tuple or list
        A tuple or list of floats indicating the degree by which the gaps in
        the environment of ascending, maximum, and descending sonority should be
        decreased or increased.

    factor : float
        A scaling factor by which the specific positions of initial and final
        should be increased and decreased.

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
    [1.5600000000000001, 1.0, 1.2, 0.69999999999999996]
    >>> prosodic_weights(prostring,scale=(4,1,2),factor=0.5)
    [6.0, 1, 4, 0.5]
    
    See also
    --------
    prosodic_string

    Todo
    ----
    Change the documentation for the new scaling mode!
    """
    
    if scale:
        if len(scale) == 3:
            print('warning, this mode is deprecated...')
            transform = {
                    '#' : scale[0] * ( 1 + factor ),
                    'V' : scale[1] * ( 1 + factor ),
                    'C' : scale[0], # consonant in ascending position
                    'c' : scale[2], # consonant in descending position
                    'v' : scale[1], # 
                    '<' : scale[1] * ( 1 - factor ), # vowel preceded by a vowel
                    '$' : scale[2] * ( 1 - factor ), # consonant in end position
                    '>' : scale[1] * ( 1 - factor ), # vowel in final position
                    'T' : scale[1],
                    '_' : 0.0
                    }
        # user-defined scale
        elif len(scale) == 10:
            transform = dict([(i[0],i[1]) for i in zip('#VCcv<$>T_',scale)])

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
                '_' : 0.0
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
                '_' : 0.0
                }
    
    out = [transform[i] for i in prostring]

    return out

def class2tokens(
        tokens,
        classes,
        gap_char = '-'
        ):
    """
    Turn aligned sound-class sequences into an aligned sequences of IPA tokens.

    Parameters
    ----------

    tokens : list
        The list of tokens corresponding to the unaligned IPA string.

    classes : string or list
        The aligned class string.

    gap_char : string
        The character which indicates gaps in the aligned class string
        (defaults to "-").
    
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

    out = []

    # copy ipa_tks in order to prevent its original form
    #_tokens = list(tokens[:])
    out = [t for t in tokens]
    for i in range(len(classes)):
        if classes[i] in ['-','X']:
            out.insert(i,'-')

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

def check_tokens(tokens):
    """
    Function checks whether tokens are given in a consistent input format.
    """
    
    errors = []

    for i,token in enumerate(tokens):
        
        # check for conversion within the articulation-model
        try:
            art.converter[token]
        except:
            try:
                art.converter[token[0]]
            except:
                errors.append((i,token))
    
    return errors

def ngrams(sequence):
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
    >>> ngrams('abcde')
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
    Convert sequence in ipa.sampa-format to ipa.unicode.
    """

    result = ''
    tokens = reXS.findall(seq)
    for tok, err in tokens:
        assert not err and tokens
        result += xs[tok]

    return result
