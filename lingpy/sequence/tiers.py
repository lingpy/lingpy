# *-* coding: utf-8 *-*
"""
Module provides tools to handle transcriptions as multi-tiered sequences.
"""
from __future__ import print_function, division, unicode_literals

from lingpy.sequence import sound_classes as sounds
from lingpy.settings import rcParams


def sound_type(sound):
    """Shortcut to determine basic sound type (C, V, or T).
    """
    return sounds.token2class(sound, 'cv')

def is_sound(sound, what):
    """
    Check whether a sound is a vowel or not.
    """
    return sound_type(sound)  == {
            'vowel': 'V', 'consonant': 'C', 'tone': 'T'}.get(what)

def is_consonant(sound):
    return is_sound(sound, 'consonant')

def is_tone(sound):
    return is_sound(sound, 'tone')

def is_vowel(sound):
    return is_sound(sound, 'vowel')

def is_stressed(sound):
    """
    Quick check for stress.
    """
    if set(rcParams['stress']).intersection(set(sound)):
        return True
    return False

def remove_stress(sound):
    return ''.join([x for x in sound if x not in rcParams['stress']])

def get_stress(sound):
    return ''.join([x for x in sound if x in rcParams['stress']])

def cvcv(sequence, **keywords):
    """
    Create a CV-template representation out of a sound sequence.
    
    """

    if not isinstance(sequence, (list, tuple)):
        sequence = sounds.ipa2tokens(sequence, **keywords)
    
    # store data in different tiers
    base, tone, stress = [], '.', []
    cv = [sound_type(sound) for sound in sequence]
    template = []
    
    # check for first instance
    if cv[0] != 'C':
        base += ['Ø.']
        template += ['Ø']

    # define preceding
    preceding = ''
    for sound_, type_ in zip(sequence, cv):
        sound = remove_stress(sound_)
        if type_ == 'C':
            if preceding == 'C':
                base[-1] += 'Ø'
                template[-1] += 'Ø'
            stress += [get_stress(sound_)]
            base += [sound+'.']
            template += ['C']
        elif type_ == 'V':
            if preceding == 'V':
                base += ['Ø.']
                template += ['Ø']
            base[-1] += sound
            template[-1] += 'V'
        elif type_ == 'T':
            type_ = preceding # nasty hack
            tone = sound

        preceding = type_


    if preceding == 'C':
        base[-1] += 'Ø'
        template[-1] += 'Ø'
    
    return base, template, stress, [tone for t in base]

if __name__ == '__main__':
    for x in ['ˈmattis', 'herbst', 'liaŋ⁵', 'hao²¹⁴', 'maˈtilde', 'geˈburt', 'ˈbraten']:
        c1, c2, c3, c4 = cvcv(x, merge_vowels=False, merge_geminates=False)
        print('\t'.join(c1))
        print('\t'.join(c2))
        print('\t'.join(c3))
        print('\t'.join(c4))
        print('')

    from lingpy.tests.util import test_data
    from lingpy import *
    from collections import defaultdict
    patterns = defaultdict(lambda : defaultdict(list))
    wl = Wordlist(test_data('KSL.qlc'))
    for idx, ipa, lang, concept in iter_rows(wl, 'ipa', 'doculect', 'concept'):
        c1, c2, c3, c4 = cvcv(ipa, merge_vowels=False, merge_geminates=False,
                expand_nasals=True)
        print('\t'.join(c1))
        patterns[lang][' '.join(c2)] += [(ipa, concept)]
    lang = 'Navajo'
    for pattern in sorted(patterns[lang]):
        print('{0:20}'.format(pattern), '\t',
                '{0:4}'.format(len(patterns[lang][pattern])), 
                ' '.join(['{0} ({1})'.format(a, b) for a, b in
                    patterns[lang][pattern][:3]]))
    print(len(patterns[lang]))

