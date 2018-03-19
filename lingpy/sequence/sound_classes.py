# *-* coding: utf-8 *-*
"""
Module provides various methods for the handling of sound classes.
"""
from __future__ import print_function, division, unicode_literals

import re
from six import text_type
import unicodedata
from collections import defaultdict

from clldutils.text import split_text, strip_brackets, strip_chars

from lingpy import log
from lingpy.util import setdefaults
from lingpy.settings import rcParams
from lingpy.data.ipa.sampa import reXS, xs


def ipa2tokens(istring, **keywords):
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

    expand_nasals : bool (default=False)

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
        breaks=rcParams['breaks'],
        combiners=rcParams['combiners'],
        diacritics=rcParams['diacritics'],
        expand_nasals=False,
        merge_geminates=True,
        merge_vowels=rcParams['merge_vowels'],
        semi_diacritics='',
        stress=rcParams['stress'],
        tones=rcParams['tones'],
        vowels=rcParams['vowels'],
        clean_sequence=False  # add this later, not today XXX
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
    vowel = False  # no vowel
    tone = False  # no tone
    merge = False  # no merge command
    start = True  # start of unit
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
        elif char in semi_diacritics and not start and not vowel \
                and not tone and out[-1] not in nogos:
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
            out += [char]
            start = False
            tone = False

    if nasal:
        out += [rcParams['nasal_placeholder']]

    if kw['merge_geminates']:
        new_out = [out[0]]
        for i in range(len(out) - 1):
            outA = out[i]
            outB = out[i + 1]
            if outA == outB:
                new_out[-1] += outB
            else:
                new_out += [outB]
        return new_out

    return out


def syllabify(seq, output="flat", **keywords):
    """
    Carry out a simple syllabification of a sequence, using sonority as a proxy.

    Parameters
    ----------
    output: {"flat", "breakpoints", "nested"} (default="flat")
        Define how to output the syllabification. Select between:
        * "flat": A syllable separator is introduced to mark the syllable boundaries
        * "breakpoins": A tuple consisting of indices that slice the original sequence into syllables is returned.
        * "nested": A nested list reflecting the syllable structure is returned.

    sep : str (default="◦")
        Select your preferred syllable separator.

    Notes
    -----

    When analyzing the sequence, we start a new syllable in all cases where we
    reach a deepest point in the sonority hierarchy of the sonority profile of
    the sequence. When passing an aligned string to this function, the gaps
    will be ignored when computing boundaries, but later on re-introduced, if
    the alignment is passed in segmented form.

    Returns
    -------

    syllable : list
        Either a flat list containing a morpheme separator, or a nested list,
        reflecting the syllable structure, or a list of tuples containing the
        indices indicating where the input sequence should be sliced in order
        to split it into syllables.

    """
    if output not in ['flat', 'breakpoints', 'nested']:
        raise ValueError(
            "The output «{0}» you specified is not available.".format(output))

    kw = {
        "sep": rcParams['morpheme_separator'],
        "gap": rcParams['gap_symbol'],
        "model": "art",
        "stress": rcParams['stress'],
        "diacritics": rcParams['diacritics'],
        "cldf": False
    }
    kw.update(keywords)

    # we assume we are dealing with tokens if the syllable is a list.
    if isinstance(seq, (list, tuple)):
        listed_seq = [s for s in seq]
    elif ' ' in seq:
        listed_seq = seq.split(' ')
    else:
        listed_seq = ipa2tokens(seq)

    # check whether our sequence is an alignment
    alm = False
    if kw['gap'] in listed_seq:
        alm = [x for x in listed_seq]
        listed_seq = [s for s in listed_seq if s != kw['gap']]

    # get the profile for the sequence
    profile = [0] + [int(i) for i in tokens2class(listed_seq, kw['model'], cldf=kw['cldf'],
        stress=kw['stress'], diacritics=kw['diacritics'])] + [0]

    new_syl = False
    breaks = []
    for i in range(1, len(profile) - 1):
        # get the pro-tokens
        p1, p2, p3 = profile[i - 1], profile[i], profile[i + 1]

        # get the char
        char = listed_seq[i - 1]

        # simple rule: we start a new syllable, if p2 is smaller or equal to p1 and p3
        # is larger than p2
        if p1 >= p2 < p3:
            if p3 == 8 or p3 == 9:
                new_syl = False
            # don't break if we are in the initial and no vowel followed
            # can be expanded to general "vowel needs to follow"-rule
            elif p1 != 7 and p2 != 7 and i == 2:
                new_syl = False

            # don't break if we are in the end of the word
            elif i == len(profile) - 3 and p3 != 7:
                new_syl = False
            else:
                new_syl = True
        # break always if there's a tone
        if p1 == 8:
            new_syl = True

        # get the char before, after
        if new_syl:
            # control for break chars which are already there
            if p1 == 9:
                breaks += [char]
            else:
                breaks += [kw['sep'], char]
            new_syl = False
        else:
            breaks += [char]

    # if we detected an alignment character in the string, we need to reparse
    # the data
    if alm:
        out = []
        idxA, idxB = 0, 0

        while idxB < len(alm) and idxA < len(breaks):
            charA = breaks[idxA]
            charB = alm[idxB]

            if charA == charB:
                idxA += 1
                idxB += 1
                out += [charA]

            elif charB == kw['gap']:
                idxB += 1
                out += [charB]

            elif charA == kw['sep']:
                idxA += 1
                out += [charA]
    else:
        out = breaks
        alm = listed_seq

    if output in ['breakpoints', 'nested']:
        bpoints = _get_breakpoints(out, kw['sep'])
        if output == 'nested':
            return [alm[x:y] for (x, y) in bpoints]
        return bpoints
    return out


def _get_breakpoints(seq, sep):
    """
    Helper function determines the points where to split a sequence.
    """
    bpoints, start, count = [], 0, 0
    for i, char in enumerate(seq):
        if char == sep:
            bpoints += [(start, count)]
            start = count
        else:
            count += 1
    return bpoints + [(start, i + 1)]


def tokens2morphemes(tokens, **keywords):
    """
    Split a string into morphemes if it contains separators.

    Notes
    -----
    Function splits a list of tokens into subsequent lists of morphemes if the list
    contains morpheme separators. If no separators are found, but tonemarkers,
    it will still split the string according to the tones. If you want to avoid
    this behavior, set the keyword **split_on_tones** to False.

    Parameters
    ----------
    sep : str (default="◦")
        Select your morpheme separator.
    word_sep: str (default="_")
        Select your word separator.

    Returns
    -------
    morphemes : list
        A nested list of the original segments split into morphemes.
    """

    if not isinstance(tokens, (list, tuple)):
        raise ValueError("The sequence needs to be a list or a tuple.")

    kw = {
        "sep": rcParams['morpheme_separator'],
        "word_sep": rcParams['word_separator'],
        "word_seps": rcParams['word_separators'],
        "seps": rcParams['morpheme_separators'],
        "split_on_tones": True,
        "tone": "T",
        "cldf": False
    }
    kw.update(keywords)
    if not kw['split_on_tones']: kw['tone'] = ''

    # check for other hints than the clean separators in the data
    new_tokens = [t for t in tokens]
    if not kw['sep'] in tokens and not kw['word_sep'] in tokens:
        class_string = tokens2class(tokens, 'cv', cldf=kw['cldf'])
        if kw['tone'] in class_string \
                and '+' not in class_string and '_' not in class_string:
            new_tokens = []
            for i, token in enumerate(tokens):
                if class_string[i] == kw['tone'] and i != len(class_string) - 1:
                    new_tokens += [token, kw['sep']]
                else:
                    new_tokens += [token]
    out = [[]]
    for i, token in enumerate(new_tokens):
        if token not in kw['sep'] + kw['word_sep'] + kw['word_seps'] + kw['seps']:
            out[-1] += [token]
        else:
            out += [[]]
    # check for bad examples
    if ['' for x in out if not x]:
        raise ValueError("[!] Your data contains empty morpheme segments.")

    return out


def _split_syllables(syllables, tokens):
    """
    Split the output of the syllabify method into subsets.

    Notes
    -----
    This is a simple helper function to deal with syllabified content.
    """

    out = [[]]
    idx = 0
    for s in syllables:
        if s != rcParams['morpheme_separator'] and s not in '#_':
            out[-1] += [(s, tokens[idx])]
            idx += 1
        else:
            # reconsider deleting these lines, since they
            # may well confuse the algorithms and we should
            # better restrict all actions to but one syllable separator
            if s in '#_':
                idx += 1
            out += [[]]

    return out


def _pprint_ono(ono):
    """
    Helper function for string output of ono-parsed words.
    """
    if rcParams['morpheme_separator'] in ono:
        return ' '.join([x[0] if not isinstance(x, text_type) else
            x for x in ono])
    out = []
    for k in ono:
        out.append(k[0] or '-')
        if k[1] == 'c':
            out.append(rcParams['morpheme_separator'])

    return ' '.join(out)


def ono_parse(word, output='', **keywords):
    """
    Carry out a rough onset-nucleus-offset parse of a word in IPA.

    Notes
    -----
    Method is an approximation and not supposed to do without flaws. It is,
    however, rather helpful in most instances. It defines a so far simple model
    in which 7 different contexts for each word are distinguished:

    * "#": onset cluster in a word's initial
    * "C": onset cluster in a word's non-initial
    * "V": nucleus vowel in a word's initial syllable
    * "v": nucleus vowel in a word's non-initial and non-final syllable
    * ">": nucleus vowel in a word's final syllable
    * "c": offset cluster in a word's non-final syllable
    * "$": offset cluster in a word's final syllable


    """
    kw = {
        "sep": rcParams['morpheme_separator'],
        "gap": rcParams['gap_symbol'],
        "model": "art",
        "stress": rcParams['stress'],
        "diacritics": rcParams['diacritics'],
        "cldf": False
    }
    kw.update(keywords)
    if isinstance(word, text_type):
        tokens = ipa2tokens(word, **kw)
    else:
        tokens = [x for x in word]
    syllabified = syllabify(tokens, **kw)
    prostring = prosodic_string(tokens, _output='CcV')
    syllables = _split_syllables(syllabified, prostring)

    out = []
    for i, syl in enumerate(syllables):

        # we take lists to restore internal tokenization
        parse = [[], [], []]
        env = 'C'
        for s in syl:
            if s[1] == 'C' and env == 'C':
                parse[0] += [s[0]]
            elif s[1] == 'C' and env != 'C':
                parse[2] += [s[0]]
            elif s[1] == 'V':
                parse[1] += [s[0]]
                env = 'V'
            elif s[1] == 'c' and env == 'C':
                parse[0] += [s[0]]

            elif s[1] == 'c':
                parse[2] += [s[0]]
                env = 'c'

        # correct parse by type
        if i == 0:
            parse[0] = (parse[0], '#')
            parse[1] = (parse[1], 'V')

            if len(syllables) != 1:
                parse[2] = (parse[2], 'c')
            else:
                parse[2] = (parse[2], '$')

        elif i == len(syllables) - 1:
            parse[0] = (parse[0], 'C')
            parse[1] = (parse[1], '>')
            parse[2] = (parse[2], '$')

        else:
            parse[0] = (parse[0], 'C')
            parse[1] = (parse[1], 'v')
            parse[2] = (parse[2], 'c')

        # linearize parse [XXX bad solution but too lazy to correct it in this
        # stage @lingulist]
        for tokens, prostring in parse:
            if tokens:
                for token in tokens:
                    out += [(token, prostring)]
            else:
                out += [('-', prostring)]

        if i < len(syllables) - 1:
            out += [rcParams['morpheme_separator']]

    if output == 'pprint':
        return _pprint_ono(out)
    elif output == 'prostring':
        return ''.join([p[1] for p in out if len(p) == 2 and p[0] != '-'])

    return out


def asjp2tokens(seq, merge_vowels=True):
    tokens = ' '.join(
        ipa2tokens(
            seq,
            diacritics='*$~"',
            vowels='aeiouE3',
            tones='',
            combiners='',
            merge_vowels=merge_vowels
        )
    )
    tokens = re.sub(r'([^ ]) ([^ ])~', r'\1\2~', tokens)
    tokens = re.sub(r'([^ ]) ([^ ]) ([^ ])\$', r'\1\2\3$', tokens)
    return tokens.split(' ')


def token2class(token, model, stress=None, diacritics=None, cldf=None):
    """
    Convert a single token into a sound-class.

    tokens : str
        A token (phonetic segment).

    model : :py:class:`~lingpy.data.model.Model`
        A :py:class:`~lingpy.data.model.Model` object.

    stress : str (default=rcParams['stress'])
        A string containing the stress symbols used in the analysis. Defaults
        to the stress as defined in ~lingpy.settings.rcParams.

    diacritics : str (default=rcParams['diacritics'])
        A string containing diacritic symbols used in the analysis. Defaults to
        the diacritic symbolds defined in ~lingpy.settings.rcParams.

    cldf : bool (default=False)
        If set to True, this will allow for a specific treatment of phonetic
        symbols which cannot be completely resolved (e.g., laryngeal h₂ in
        Indo-European). Following the `CLDF <http://cldf.clld.org>`_
        specifications (in particular the specifications for writing
        transcriptions in segmented strings, as employed by the `CLTS
        <http://calc.digling.org/clts/>`_ initiative), in cases of insecurity
        of pronunciation, users can adopt a ```source/target``` style, where
        the source is the symbol used, e.g., in a reconstruction system, and
        the target is a proposed phonetic interpretation. This practice is also
        accepted by the `EDICTOR <http://edictor.digling.org>`_ tool.

    Returns
    -------

    sound_class : str
        A sound-class representation of the phonetic segment. If the segment
        cannot be resolved, the respective string will be rendered as "0"
        (zero).

    See also
    --------
    ipa2tokens
    class2tokens
    token2class

    """
    # check basic parameters
    stress = rcParams['stress'] or stress
    diacritics = rcParams['diacritics'] or diacritics

    # change token if cldf is selected
    if cldf:
        token = token.split('/')[1] or '?' if '/' in token else token

    # check whether model is passed as real model or as string
    if str(model) == model:
        model = rcParams[model]

    try:
        return model[token]
    except KeyError:
        try:
            return model[token[0]]
        except IndexError:
            return '0'
        except KeyError:
            # check for stressed syllables
            if len(token) > 0:
                if token[0] in stress and len(token) > 1:
                    try:
                        return model[token[1:]]
                    except KeyError:
                        try:
                            return model[token[1]]
                        except KeyError:
                            # new character for missing data and spurious items
                            return '0'
                elif token[0] in diacritics:
                    if len(token) > 1:
                        try:
                            return model[token[1:]]
                        except KeyError:
                            try:
                                return model[token[1]]
                            except KeyError:
                                return '0'
                    else:
                        return '0'
                else:
                    # new character for missing data and spurious items
                    return '0'
            else:
                return '0'


def tokens2class(tokens, model, stress=None, diacritics=None, cldf=False):
    """
    Convert tokenized IPA strings into their respective class strings.

    Parameters
    ----------

    tokens : list
        A list of tokens as they are returned from :py:func:`ipa2tokens`.

    model : :py:class:`~lingpy.data.model.Model`
        A :py:class:`~lingpy.data.model.Model` object.

    stress : str (default=rcParams['stress'])
        A string containing the stress symbols used in the analysis. Defaults
        to the stress as defined in ~lingpy.settings.rcParams.

    diacritics : str (default=rcParams['diacritics'])
        A string containing diacritic symbols used in the analysis. Defaults to
        the diacritic symbolds defined in ~lingpy.settings.rcParams.

    cldf : bool (default=False)
        If set to True, this will allow for a specific treatment of phonetic
        symbols which cannot be completely resolved (e.g., laryngeal h₂ in
        Indo-European). Following the `CLDF <http://cldf.clld.org>`_ specifications (in particular the
        specifications for writing transcriptions in segmented strings, as
        employed by the `CLTS <http://calc.digling.org/clts/>`_ initiative), in
        cases of insecurity of pronunciation, users can adopt a
        ```source/target``` style, where the source is the symbol used, e.g.,
        in a reconstruction system, and the target is a proposed phonetic
        interpretation. This practice is also accepted by the `EDICTOR
        <http://edictor.digling.org>`_ tool.

    Returns
    -------

    classes : list
        A sound-class representation of the tokenized IPA string in form of a
        list. If sound classes cannot be resolved, the respective string will
        be rendered as "0" (zero).

    Notes
    -----
    The function ~lingpy.sequence.sound_classes.token2class returns a "0"
    (zero) if the sound is not recognized by LingPy's sound class models. While
    an unknown sound in a longer sequence is no problem for alignment
    algorithms, we have some unwanted and often even unforeseeable behavior,
    if the sequence is completely unknown. For this reason, this function
    raises a ValueError, if a resulting sequence only contains unknown sounds.

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
    token2class

    """
    # raise value error if input is not an iterable (tuple or list)
    if not isinstance(tokens, (tuple, list)):
        raise ValueError("[!] Need tuple or list as input.")

    stress=rcParams['stress']
    diacritics=rcParams['diacritics']

    out = []
    for token in tokens:
        out.append(token2class(token, model, stress=stress,
            diacritics=diacritics, cldf=cldf))
    if out.count('0') == len(out):
        raise ValueError("[!] your sequence contains only unknown characters")
    return out


def prosodic_string(string, _output=True, **keywords):
    """
    Create a prosodic string of the sonority profile of a sequence.

    Parameters
    ----------

    seq : list
        A list of integers indicating the sonority of the tokens of the
        underlying sequence.

    stress : str (default=rcParams['stress'])
        A string containing the stress symbols used in the analysis. Defaults
        to the stress as defined in ~lingpy.settings.rcParams.

    diacritics : str (default=rcParams['diacritics'])
        A string containing diacritic symbols used in the analysis. Defaults to
        the diacritic symbolds defined in ~lingpy.settings.rcParams.

    cldf : bool (default=False)
        If set to True, this will allow for a specific treatment of phonetic
        symbols which cannot be completely resolved (e.g., laryngeal h₂ in
        Indo-European). Following the `CLDF <http://cldf.clld.org>`_
        specifications (in particular the specifications for writing
        transcriptions in segmented strings, as employed by the `CLTS
        <http://calc.digling.org/clts/>`_ initiative), in cases of insecurity
        of pronunciation, users can adopt a ```source/target``` style, where
        the source is the symbol used, e.g., in a reconstruction system, and
        the target is a proposed phonetic interpretation. This practice is also
        accepted by the `EDICTOR <http://edictor.digling.org>`_ tool.

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
    * ``T``: tone
    * ``_``: word break

    Examples
    --------
    >>> prosodic_string(ipa2tokens('t͡sɔyɡə')
    'AXBZ'

    """
    defaults = dict(stress=rcParams['stress'], cldf=False,
        diacritics=rcParams['diacritics'])
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # check for empty string passed
    if not string:
        return ''

    # check for the right string
    if not all(isinstance(x, int) for x in string):
        # get the sonority profile
        sstring = [9] + \
                  [int(t) for t in tokens2class(string, rcParams['art'],
                      stress=keywords['stress'],
                      diacritics=keywords['diacritics'], cldf=keywords['cldf'])] + \
                  [9]
    else:
        sstring = [9] + string + [9]

    # check for multiple strings in string
    if 9 in sstring[1:-1]:
        # break the string into pieces
        nstrings = [[]]
        for i in sstring[1:-1]:
            if i != 9:
                nstrings[-1] += [i]
            else:
                nstrings += [[]]

        # return the prostrings of the pieces recursively, note that the
        # additional check whether x is True is necessitated by the fact that
        # often errors occur in the coding, i.e. strings are given
        return '_'.join(prosodic_string(x, _output) for x in nstrings)

    # create the output values
    pstring = ''
    first = True  # stores whether first syllable is currently being processed

    # start iterating over relevant parts of the string
    for i in range(1, len(sstring) - 1):
        # get all three values
        a, b, c = sstring[i - 1], sstring[i], sstring[i + 1]

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
                        pstring = pstring[:-1] + pstring[-1].replace('L', 'M') + 'B'
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
            raise ValueError(
                "Conversion to prosodic string failed due to a condition which was not "
                "defined in the convertion, for details compare the numerical string "
                "{0} with the profile string {1}".format(sstring, pstring))

    if _output == 'cv':
        conv = {
            "A": "C",
            "B": "C",
            "C": "C",
            "M": "C",
            "L": "C",
            "N": "C",
            "X": "V",
            "Y": "V",
            "Z": "V",
            "T": "T",
            "_": "_",
        }
        return ''.join([conv[x] for x in pstring])

    elif _output == 'CcV':
        conv = {
            "A": "C",
            "B": "C",
            "C": "C",
            "M": "c",
            "L": "c",
            "N": "c",
            "X": "V",
            "Y": "V",
            "Z": "v",
            "T": "T",
            "_": "_",
        }
        return ''.join([conv[x] for x in pstring])

    elif _output:
        return pstring

    else:
        conv = {
            "A": "#",
            "B": "C",
            "C": "C",
            "M": "c",
            "L": "c",
            "N": "$",
            "X": "V",
            "Y": "v",
            "Z": ">"
        }
        return ''.join([conv[x] for x in pstring])


def prosodic_weights(prostring, _transform={}):
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
            '#': 1.6,
            'V': 3.0,
            'c': 1.1,
            'v': 3.0,  # raise the cost for the gapping of vowels
            '<': 0.8,
            '$': 0.5,
            '>': 0.7,

            # new values for alternative prostrings
            'A': 1.6,  # initial
            'B': 1.3,  # syllable-initial
            'C': 1.2,  # ascending

            'L': 1.1,  # descending
            'M': 1.1,  # syllable-descending
            'N': 0.5,  # final

            'X': 3.0,  # vowel in initial syllable
            'Y': 3.0,  # vowel in non-final syllable
            'Z': 0.7,  # vowel in final syllable
            'T': 1.0,  # Tone
            '_': 0.0  # break character
        }
    # default scale for other languages
    else:
        transform = {
            '#': 2.0,
            'V': 1.5,
            'c': 1.1,
            'v': 1.3,
            '<': 0.8,
            '$': 0.8,
            '>': 0.7,

            # new values for alternative prostrings
            'A': 2.0,  # initial
            'B': 1.75,  # syllable-initial
            'C': 1.5,  # ascending

            'L': 1.1,  # descending
            'M': 1.1,  # syllable-descending
            'N': 0.8,  # final

            'X': 1.5,  # vowel in initial syllable
            'Y': 1.3,  # vowel in non-final syllable
            'Z': 0.8,  # vowel in final syllable
            'T': 0.0,  # Tone
            '_': 0.0  # break character

        }

    return [transform[i] for i in prostring]


def class2tokens(tokens, classes, gap_char='-', local=False):
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
        for i, c in enumerate(classes):
            if c in '-X':
                out.insert(i, gap_char)
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
                out.insert(i, gap_char)
    return out


def pid(almA, almB, mode=2):
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

    zipped = zip(almA, almB)
    idn_pos = 0
    int_gps = 0
    aln_pos = 0

    for charA, charB in zipped:
        tmp = [charA, charB].count('-')
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
            log.warn('Zero Division Error in {0} and {1}'.format(almA, almB))
            return 0

    elif mode == 1:
        try:
            return idn_pos / aln_pos
        except ZeroDivisionError:
            log.warn('Zero Division Error in {0} and {1}'.format(almA, almB))
            return 0

    elif mode == 3:
        srt_seq = min(
            len([i for i in almA if i != '-']), len([i for i in almB if i != '-']))
        try:
            return idn_pos / srt_seq
        except ZeroDivisionError:
            log.warn('Zero Division Error in {0} and {1}'.format(almA, almB))
            return 0

    elif mode == 4:
        srt_seq = min(
            len(''.join([i[0] for i in almA]).strip('-')),
            len(''.join([i[0] for i in almB]).strip('-')))
        try:
            return idn_pos / srt_seq
        except ZeroDivisionError:
            log.warn('Zero Division Error in {0} and {1}'.format(almA, almB))
            return 0

    elif mode == 5:
        return idn_pos / len(almA)

def check_tokens(tokens, **keywords):
    """
    Function checks whether tokens are given in a consistent input format.
    """
    setdefaults(keywords, stress=rcParams['stress'],
            diacritics=rcParams['diacritics'], cldf=False)
    errors = []
    for i, token in enumerate(tokens):
        # check for conversion within the articulation-model
        cls = token2class(token, rcParams['art'], stress=keywords['stress'],
                cldf=keywords['cldf'], diacritics=keywords['diacritics'])
        if cls == '0':
            errors.append((i, token))

    return errors

def get_all_ngrams(sequence, sort=False):
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
        for j in range(1, len(new_sequence)):
            out += [new_sequence[:j]]
            out += [new_sequence[j:]]

        # increment i and decrement l
        i += 1
        l -= 1

    sort = sort or list

    return sort(out)


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
        try:
            assert not err and tokens
        except AssertionError:
            raise AssertionError('{0} + {1}'.format(err, tok))
        result += xs[tok]

    return result


def _seq_as_list(sequence):
    if ' ' in sequence and isinstance(sequence, str):
        return sequence.split(' ')
    return list(sequence)


def bigrams(sequence):
    """
    Convert a given sequence into a sequence of bigrams.
    """
    seq = _seq_as_list(sequence)
    return list(zip(['#'] + seq, seq + ['$']))


def trigrams(sequence):
    """
    Convert a given sequence into a sequence of trigrams.
    """
    seq = _seq_as_list(sequence)
    return list(zip(['#', '#'] + seq, ['#'] + seq + ['$'], seq + ['$', '$']))


def fourgrams(sequence):
    """
    Convert a given sequence into a sequence of trigrams.
    """
    seq = _seq_as_list(sequence)
    return list(
        zip(
            ['#', '#', '#'] + seq,
            ['#', '#'] + seq + ['$'],
            ['#'] + seq + ['$', '$'],
            seq + ['$', '$', '$']
        )
    )


def get_n_ngrams(sequence, ngram=4):
    """
    convert a given sequence into a sequence of ngrams.
    """
    seq = _seq_as_list(sequence)
    tobezipped = []
    for i in range(ngram):
        prefix = (ngram - i - 1) * ['#']
        postfix = i * ['$']
        tobezipped += [tuple(prefix + seq + postfix)]

    return list(zip(*tobezipped))[ngram - 1:]


def pgrams(sequence, **keywords):
    """
    Convert a given sequence into bigrams consisting of prosodic string symbols and the
    tokens of the original sequence.
    """
    if isinstance(sequence, str) and ' ' not in sequence:
        seq = ipa2tokens(sequence)
    else:
        seq = _seq_as_list(sequence)
    return list(zip(seq, prosodic_string(seq, **keywords)))

def _get_brackets(brackets):

    out = defaultdict(str)
    for b in brackets:
        out[b] = unicodedata.lookup(unicodedata.name(b).replace('LEFT', 'RIGHT'))
        if b == out[b]:
            log.warn('lingpy.sequence.sound_classes.get_brackets' + \
                    'Item «{0}» does not have a counterpart!'.format(b))
    return out

def clean_string(
        sequence, semi_diacritics='hsʃ̢ɕʂʐʑʒw', merge_vowels=False,
        segmentized=False, rules=None, ignore_brackets=True, brackets=None,
        split_entries=True, splitters='/,;~', preparse=None,
        merge_geminates=True, normalization_form="NFC"):
    """
    Function exhaustively checks how well a sequence is understood by \
            LingPy.

    Parameters
    ----------
    semi_diacritics : str
        Indicate characters which can occur both as "diacritics" (second part
        in a sound) or alone.
    merge_vowels : bool (default=True)
        Indicate whether consecutive vowels should be merged.
    segmentized : False
        Indicate whether the input string is already segmentized or not. If set
        to True, items in brackets can no longer be ignored.
    rules : dict
        Replacement rules to be applied to a segmentized string.
    ignore_brackets : bool
        If set to True, ignore all content within a given bracket.
    brackets : dict
        A dictionary with opening brackets as key and closing brackets as
        values. Defaults to a pre-defined set of frequently occurring brackets.
    split_entries : bool (default=True)
        Indicate whether multiple entries (with a comma etc.) should be split
        into separate entries.
    splitters : str
        The characters which force the automatic splitting of an entry.
    prepares : list
        List of tuples, giving simple replacement patterns (source and target),
        which are applied before every processing starts.

    Returns
    -------
    cleaned_strings : list
        A list of cleaned strings which are segmented by space characters. If
        splitters are encountered, indicating that the entry contains two
        variants, the list will contain one for each element in a separate
        entry. If there are no splitters, the list has only size one.
    """
    sequence = unicodedata.normalize(normalization_form, sequence)
    rules = rules or {}
    preparse = preparse or []

    # replace white space if not indicated otherwise
    if segmentized:
        segment_list = [sequence.split(' ') if not isinstance(sequence, (list,
            tuple)) else sequence]
    else:
        for s, t in preparse:
            sequence = sequence.replace(s, t)
        segment_list = []
        if ignore_brackets:
            new_sequence = strip_brackets(sequence, brackets=brackets)
        else:
            new_sequence = sequence

        # splitting needs to be done afterwards
        if split_entries:
            new_sequences = split_text(new_sequence, splitters,
                    brackets='' if not ignore_brackets else brackets)
        else:
            new_sequences = [new_sequence]

        for new_sequence in new_sequences:
            segments = ipa2tokens(
                    re.sub(r'\s+', '_', new_sequence.strip()),
                    semi_diacritics=semi_diacritics,
                    merge_vowels=merge_vowels,
                    merge_geminates=merge_geminates)
            segment_list += [segments]
    out = []
    for segments in segment_list:
        segments = [rules.get(s, s) for s in segments]
        out += [' '.join(segments)]
    return out

def codepoint(s):
    "Return unicode codepoint(s) for a character set."
    return ' '.join(['U+'+('000'+hex(ord(x))[2:])[-4:] for x in s])


