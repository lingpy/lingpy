# *-* coding: utf-8 *-*
"""
This module provides miscellaneous functions which are mostly used internally.
"""
from __future__ import division,print_function
from re import sub,findall
from numpy import array,sqrt,zeros
from ..data import *

def ipa2tokens(
        seq, 
        diacritics = None, 
        vowels = None,
        merge_vowels = True
        ):
    """
    Tokenize IPA-encoded strings. 
    
    Parameters
    ----------

    seq : string or unicode
        The input sequence that shall be tokenized.
    
    diacritics : unicode
        A string containing all diacritics which shall be considered in the
        respective analysis. When set to *None*, the default diacritic string
        will be used.
    
    vowels : unicode
        A string containing all vowel symbols which shall be considered in the
        respective analysis. When set to *None*, the default vowel string will
        be used.
    
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
    [u't\u0361s', u'\u0254y', u'\u0261', u'\u0259']
    >>> for t in ipa2tokens(myseq): print t
    t͡s	
    ɔy	
    ɡ	
    ə
    
    See also
    --------
    tokens2class
    class2tokens

    """
    
    try:
        seq = unicode(seq,'utf-8')
    except:
        pass

    # if vowels and diacritics are not defined, take the default values
    if not diacritics:
        diacritics = ipa_diacritics
    if not vowels:
        vowels = ipa_vowels
    
    # merge vowels and diacritics for later use
    dv = diacritics + vowels
    tones = unicode('⁰¹²³⁴⁵⁶⁷⁸⁹⁻₀₁₂₃₄₅₆','utf-8')
    
    # replace all multiple dots by just one dot (needed for local alignment
    # analyses) in order to keep trace of dotted sequences
    seq = sub(r'\.+',r'.',seq)

    # replace double consonants by the ipa-character for long consonants
    seq = sub(r'([^'+dv+tones+r'])\1+',r'\1'+unicode('ː','utf-8'),seq)

    # matches for regular expressions using the findall-function
    # m_1 matches all affricates and following diacritics
    m_1 = r'([^'+dv+r']['+u'\u0361\u035c'+r'][^'+dv+r']['+diacritics+r']*)'
    
    # m_2 matches all vowels and following diacritics. If merge_vowels is
    # set to False, the vowels will be separated
    if merge_vowels == False:
        m_2 = r'(['+vowels+r']['+diacritics+r']*)'
    else:
        m_2 = r'(['+vowels+r']+['+dv+r']*)'
    
    # matches for consonants and diacritics
    m_3 = r'([^'+dv+tones+r']['+diacritics+r']*)'

    # matches for tones
    m_4 = r'(['+tones+r']+)' #['+tones+r']*)'

    # carry out the search. The findall function returns touples of three
    # elements, since three conditions are being checked. 
    match_list = findall(m_1+'|'+m_2+'|'+m_3+'|'+m_4,seq)
    
    # extract the non-empty matches from the match_list. The dot-character
    # (.) as a character for syllable breaks is hereby ignored.
    tokens = []
    for m1,m2,m3,m4 in match_list:
        if '.' not in [m1,m2,m3,m4]:
            if m1 != '':
                tokens.append(m1)
            elif m2 != '':
                tokens.append(m2)
            elif m3 != '':
                tokens.append(m3)
            else:
                tokens.append(m4)
    
    # check, whether input and output are of the same length, in order to
    # make sure that the regex matched all conditions one is facing
    if len(''.join(tokens)) != len(seq.replace('.','')):
        print('_'.join(tokens),len(tokens),seq.replace('.',''),len(seq.replace('.','')))
        print("[!] The length of input and output string is different.",)
        print("... There might be a coding problem.")

    return tokens

def tokens2class(
        tokens,
        model
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
    for token in tokens:
        try:
            out.append(model.converter[token])
        except:
            try:
                out.append(model.converter[token[0]])
            except:
                # new character for missing data and spurious items
                out.append('0')

    return ''.join(out)

def prosodic_string(
        seq,
        prs_model = 4
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
    
    out = [0 for i in seq]
    
    # vowel in initial position
    if seq[0] == 7:
        out[0] = 'V'
        up = True

    # consonant in initial position
    else:
        #out[0] = 'C'
        out[0] = '#'
        up = False

    if len(seq) == 1:
        return ''.join(out)
    
    i = 1
    while i < len(seq) -1:
        
        # if the item is a missing data character, we simply assign it the same
        # value as we would with a vowel (note that python allows to
        # compare strings with integers, so we use that functionality here, ?
        # scores higher than 10)
        if seq[i] == '0':
            out[i] = 'v'
            up = True
            i += 1
        # if the item is a break
        elif seq[i] == 9:
            out[i] = '_'
            up = True
            i += 1
        # if the item is a tone
        elif seq[i] == 8:
            out[i] = 'T'
            ## change outgoing consonants in tonal strings to vowels
            #if out[i-1] == 'c':
            #    out[i-1] = '<'
            up = True
            i += 1
        # a vowel in the first position in a syllable XXX has just been added
        # and has to be checked, whether this line is useful XXX
        #elif seq[i] == 7 and 'V' not in out[:i]:
        #    out[i] = 'V'
        #    up = False
        #    i += 1
        
        # a vowel is not preceded by a vowel
        elif seq[i] == 7 and seq[i-1] != 7:
            out[i] = 'v'
            up = False
            i += 1
        # a vowel is preceded by a vowel
        elif seq[i] == 7:
            out[i] = '<'
            up = False
            i += 1
        ## XXX a consonant in syllable-initial position XXX
        #elif seq[i] < 7 and out[i-1] == '#':
        #    out[i] = '#'
        #    i += 1

        # a consonant in ascending position
        elif seq[i] <= seq[i+1] and up and seq[i+1] != 8:
            out[i] = 'C'
            i += 1
        # a consonant in descending position
        elif not up and seq[i] >= seq[i+1] or seq[i+1] == 8:
            out[i] = 'c'
            i += 1
        else:
            out[i] = 'C'
            i += 1
            up = True

    # get the last position
    if seq[-1] == 7: 
        out[-1] = '>' # vowel in final position
    elif seq[-1] == 8:
        out[-1] = 'T' # tone in final position
        ## change outgoing consonants in tonal strings to vowels
        #if out[i-1] == 'c':
        #    out[i-1] = '<'
    else:
        out[-1] = '$' # consonant in final position
        #out[-1] = 'c'
    
    # XXX add-on simplifies the handling of alternative models via the use of
    # regexes XXX
    out = ''.join(out)
    if prs_model == 1:
        return out
    elif prs_model == 2:
        out = out.replace('#','C')
        return out
    elif prs_model == 3:
        v = out.find('v')
        if v:
            out = out[:v].replace('C','#')+out[v:]
        c = out.find('v')
        if c:
            out = out[:c+1].replace('v','V')+out[c+1:]
        return out
    elif prs_model == 4:
        #v = out.find('v')
        #if v:
        #    out = out[:v].replace('C','#')+out[v:]
        c = out.find('v')
        if c:
            out = out[:c+1].replace('v','V')+out[c+1:]
        return out


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
    #i = len(classes) - 1
    #while i >= 0:
    #    if classes[i] not in ['-', 'X']:
    #        out = [_tokens.pop()] + out 
    #        i-=1
    #    else:
    #        out = [gap_char] + out
    #        i-=1

    return out

def squareform(x):
    """
    A simplified version of the :py:func:`scipy.spatial.distance.squareform` \
    function.

    Parameters
    ----------

    x : :py:class:`numpy.array` or list
        The one-dimensional flat representation of a symmetrix distance matrix.

    Returns
    -------
    matrix : :py:class:`numpy.array`
        The two-dimensional redundant representation of a symmetric distance matrix.

    """

    l = len(x)

    # calculate the length of the square
    s = int(sqrt(2 * l) + 1)
    
    out = zeros((s,s))
    
    k = 0
    for i in range(s):
        for j in range(s):
            if i < j:
                out[i][j] = x[k]
                out[j][i] = x[k]
                k += 1
    return out

def loadtxt(infile):
    """
    Function imitates the :py:func:`numpy.loadtxt` function.

    Parameters
    ----------
    infile : file
        The input file from which the data is read.

    Returns
    -------
    data : list
        A list object which renders the dimensions of the input file.

    """

    f = open(infile)

    data = []
    for line in f:
        if not line.startswith('#'):
            data.append([x.strip() for x in line.strip().split('\t')])

    # check the data for consistency
    start = len(data[0])
    problems = False
    for i,x in enumerate(data):
        if len(x) != start:
            problems = True
            break

    # if there are inconsistent lines, an integrity check of the data has to be
    # carried out, ask the user whether this is wanted, or not
    if problems:
        answer = raw_input('[!] Input file contains errors. Do you want the errors to be listed separatedly (y/n)? ')
        if answer == 'y':
            for i,x in enumerate(data):
                if len(x) != start:
                    print("... {0} columns in line {1}, expected {2} ...\
                            ...".format(
                                len(x),
                                i+1,
                                start
                                )
                            )
        else:
            pass
    
    return data


class LingpyArray(object):
    """
    An extension of the numpy array object which allows the storage of lists in
    two-dimensional arrays.

    Parameters
    ----------

    input_list : list
        The list which shall be converted in an array-like object.

    """

    def __init__(self,input_list):

        
        h = len(input_list)
        w = len(input_list[0])

        self.array = zeros((h,w),dtype='int')

        self.dictionary = {}

        count = 0
        for i,line in enumerate(input_list):
            for j,itm in enumerate(line):
                self.dictionary[count] = input_list[i][j]
                self.array[i][j] = count
                count += 1

    def __getitem__(self,x):

        ind = self.array[x]

        try:
            return self.dictionary[self.array[x]]
        except:
            try:
                return [self.dictionary[i] for i in ind]
            except:
                return [[self.dictionary[i] for i in j] for j in ind]


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

