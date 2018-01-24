# *-* coding: utf-8 *-*
"""
Module provides methods for the handling of orthography profiles.
"""
from __future__ import print_function, division, unicode_literals

from six import text_type
import unicodedata
from collections import defaultdict

from clldutils.text import split_text, strip_brackets, strip_chars

from lingpy.sequence.sound_classes import codepoint, clean_string, token2class

def simple_profile(wordlist, ref='ipa', semi_diacritics='hsʃ̢ɕʂʐʑʒw', merge_vowels=False,
        brackets=None, splitters='/,;~', merge_geminates=True,
        bad_word="<???>", bad_sound="<?>", clts=None, unknown_sound="!{0}"):
    """
    Create an initial Orthography Profile using Lingpy's clean_string procedure.

    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        A wordlist from which you want to derive an initial
        orthography profile.
    ref : str (default="ipa")
        The name of the reference column in which the words are stored.
    semi_diacritics : str
        Indicate characters which can occur both as "diacritics" (second part
        in a sound) or alone.
    merge_vowels : bool (default=True)
        Indicate whether consecutive vowels should be merged.
    brackets : dict
        A dictionary with opening brackets as key and closing brackets as
        values. Defaults to a pre-defined set of frequently occurring brackets.
    splitters : str
        The characters which force the automatic splitting of an entry.
    clts : dict (default=None)
        A dictionary(like) object that converts a given source sound into a
        potential target sound, using the get()-method of the dictionary.
        Normally, we think of a CLTS instance here (that is: a cross-linguistic
        transcription system as defined in the pyclts package).
    bad_word : str (default="«???»")
        Indicate how words that could not be parsed should be handled. Note
        that both "bad_word" and "bad_sound" are format-strings, so you can add
        formatting information here.
    bad_sound : str (default="«?»")
        Indicate how sounds that could not be converted to a sound class be
        handled. Note that both "bad_word" and "bad_sound" are format-strings,
        so you can add formatting information here.
    unknown_sound : str (default="!{0}")
        If with_clts is set to True, use this string to indicate that sounds
        are classified as "unknown sound" in the CLTS framework.    
    
    Returns
    -------
    profile : generator
        A generator of tuples (three items), indicating the segment, its frequency,
        the conversion to sound classes in the Dolgopolsky sound-class model,
        and the unicode-codepoints.
    """
    clts = clts or {}
    nulls = set()
    bad_words = set()
    brackets = brackets or "([{『（₍⁽«)]}）』⁾₎"
    profile = defaultdict(int)
    words = [wordlist[idx, ref] for idx in wordlist]
    for word in words:
        if isinstance(word, list):
            word = ' '.join(word)
        cleaned_string = clean_string(word, semi_diacritics=semi_diacritics,
                merge_vowels=merge_vowels, brackets=None, ignore_brackets=False,
                split_entries=False, preparse=None, rules=None,
                merge_geminates=merge_geminates)[0]

        # retain whole word if there are splitters in the word
        if [x for x in cleaned_string if x in brackets + splitters]:
            profile[word] += 1
            bad_words.add(word)
        else:
            for segment in cleaned_string.split(' '):
                profile[segment] += 1
            for segment in [x for x in word if x not in cleaned_string]:
                profile[segment] += 1
                nulls.add(segment)

    for s, f in sorted(profile.items(), key=lambda x: x[1], reverse=True):
        sclass = token2class(s, 'dolgo')
        if s in bad_words:
            ipa = bad_word.format(s)
        elif sclass == '0' and s not in nulls:
            ipa = bad_sound.format(s)
        elif s in nulls:
            ipa = 'NULL'
        elif clts:
            sound = clts.get(s, False)
            if not sound:
                ipa = '!'+s
            else:
                ipa = text_type(sound)
        else:
            ipa = s
        yield s, ipa, text_type(f), codepoint(s)

def context_profile(wordlist, ref='ipa', col="doculect",
        semi_diacritics='hsʃ̢ɕʂʐʑʒw', merge_vowels=False, brackets=None,
        splitters='/,;~', merge_geminates=True, clts=False,
        bad_word="<???>", bad_sound="<?>", unknown_sound="!{0}", examples=2):
    """
    Create an advanced Orthography Profile with context and doculect information.

    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        A wordlist from which you want to derive an initial
        orthography profile.
    ref : str (default="ipa")
        The name of the reference column in which the words are stored.
    col : str (default="doculect")
        Indicate in which column the information on the language variety is
        stored.
    semi_diacritics : str
        Indicate characters which can occur both as "diacritics" (second part
        in a sound) or alone.
    merge_vowels : bool (default=True)
        Indicate whether consecutive vowels should be merged.
    brackets : dict
        A dictionary with opening brackets as key and closing brackets as
        values. Defaults to a pre-defined set of frequently occurring brackets.
    splitters : str
        The characters which force the automatic splitting of an entry.
    clts : dict (default=None)
        A dictionary(like) object that converts a given source sound into a
        potential target sound, using the get()-method of the dictionary.
        Normally, we think of a CLTS instance here (that is: a cross-linguistic
        transcription system as defined in the pyclts package).
    bad_word : str (default="«???»")
        Indicate how words that could not be parsed should be handled. Note
        that both "bad_word" and "bad_sound" are format-strings, so you can add
        formatting information here.
    bad_sound : str (default="«?»")
        Indicate how sounds that could not be converted to a sound class be
        handled. Note that both "bad_word" and "bad_sound" are format-strings,
        so you can add formatting information here.
    unknown_sound : str (default="!{0}")
        If with_clts is set to True, use this string to indicate that sounds
        are classified as "unknown sound" in the CLTS framework.
    examples : int(default=2)
        Indicate the number of examples that should be printed out.

    Returns
    -------
    profile : generator
        A generator of tuples (three items), indicating the segment, its frequency,
        the conversion to sound classes in the Dolgopolsky sound-class model,
        and the unicode-codepoints.
    """
    clts_ = clts or {}
    nulls = set()
    bad_words = set()
    brackets = brackets or "([{『（₍⁽«)]}）』⁾₎"
    profile = defaultdict(list)
    for idx, word, language in wordlist.iter_rows(ref, col):
        if isinstance(word, list):
            word = ' '.join(word)
        cleaned_string = clean_string(word, semi_diacritics=semi_diacritics,
                merge_vowels=merge_vowels, brackets=None, ignore_brackets=False,
                split_entries=False, preparse=None, rules=None,
                merge_geminates=merge_geminates)[0].split(' ')

        # retain whole word if there are splitters in the word
        if [x for x in cleaned_string if x in brackets + splitters]:
            profile[word] += [(language, word)]
            bad_words.add(word)
        else:
            context_pre = ['^'] + (len(cleaned_string) - 1) * ['']
            context_post = (len(cleaned_string)-1) * [''] + ['$']
            for ctxA, ctxB, segment in zip(context_pre, context_post, cleaned_string):
                profile[ctxA+segment+ctxB] += [(language, word)]
            for segment in [x for x in word if x not in ' '.join(cleaned_string)]:
                profile[segment] += [(language, word)]
                nulls.add(segment)
    
    for s in '^$':
        yield s, 'NULL', '', '', '', ''

    for idx, (s, entries) in enumerate(sorted(profile.items(), key=lambda x:
        len(x[1]), reverse=True)):
        sclass = token2class(s.strip('^$'), 'dolgo')
        words, langs = [l[1] for l in entries], [l[0] for l in entries]
        languages = ', '.join(sorted(set(langs), key=lambda x: langs.count(x),
            reverse=True))
        frequency = str(len(langs))
        codepoints = codepoint(s)
        examples_ = ', '.join(sorted(set(words), key=lambda x:
            words.count(x), reverse=True)[:examples])
        if s in bad_words:
            ipa = bad_word.format(s)
        elif sclass == '0':
            ipa = bad_sound.format(s)
        elif s in nulls:
            ipa = 'NULL'
        elif clts_:
            sound = clts_.get(s.strip('^$'), False)
            if not sound:
                ipa = '!'+s.strip('^$')
            else:
                ipa = text_type(sound)
        else:
            ipa = s.strip('^$')
        
        yield s, ipa, examples_, languages, frequency, codepoints
