# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-17 10:40
# modified : 2014-12-07 13:32
"""
Module handels all global parameters used in a LingPy session.
"""

__author__="Johann-Mattis List"
__date__="2014-12-07"

# builtin imports
from datetime import datetime,date
import os

# internal imports
from ._settings import rcParams
from .data.model import Model, load_dvt

# load diacritics, vowels, tones
diacritics, vowels, tones = load_dvt()

# these are lexstat-specific parameters, all prefixed by "lexstat"
lexstat = dict(
        lexstat_transform      = {
                    'A':'C',
                    'B':'C',
                    'C':'C',
                    'L':'c',
                    'M':'c',
                    'N':'c',
                    'X':'V', #
                    'Y':'V', #
                    'Z':'V', #
                    'T':'T', #
                    '_':'_'
                    },
        lexstat_runs           = 1000,
        lexstat_modes          = [("global",-2,0.5),("local",-1,0.5)],
        lexstat_rands          = 1000,
        lexstat_limit          = 10000,
        lexstat_scoring_method = 'shuffle',
        lexstat_ratio          = (2,1),
        lexstat_vscale         = 0.5,
        lexstat_threshold      = 0.3,
        lexstat_cluster_method = 'upgma',
        lexstat_preprocessing_method = 'sca',
        lexstat_preprocessing_threshold = 0.7
        )
rcParams.update(lexstat)

# these are alignment-specific parameters, all prefixed by "align"
alignments = dict(
        align_mode                                    = 'global',
        align_modes                                   = [
                ('global',-2,0.5),
                ('local',-1,0.5),
                ],
        align_scale                                   = 0.5,
        align_factor                                  = 0.3,
        align_gap_weight                              = 0.5,
        align_classes                                 = True,
        align_sonar                                   = True,
        align_scorer                                  = {},
        align_tree_calc                               = 'neighbor',
        align_gop                                     = -2,
        align_transform                               = {

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
            },
        align_notransform                             = {
            # new values for alternative prostrings
            'A' : 1,  # initial
            'B' : 1, # syllable-initial
            'C' : 1,  # ascending
            'L' : 1,  # descending
            'M' : 1,  # syllable-descending
            'N' : 1,  # final
            'X' : 1,  # vowel in initial syllable
            'Y' : 1,  # vowel in non-final syllable
            'Z' : 1,  # vowel in final syllable
            'T' : 1,  # Tone
            '_' : 1   # break character
            },
        align_stamp                                   = """# MSA
# dataset    : {0}
# collection : {1} 
# aligned by : LingPy-2.2 <www.lingpy.org>
# created on : {2}
# parameters : {3}
"""

        )
rcParams.update(alignments)

# dictionary stores basic parameters that are used during a LingPy session
rcParamsUpd = dict(
        verbose                    = False,
        debug                      = False,
        schema                     = 'qlc',
        asjp                       = Model('asjp'),
        sca                        = Model('sca'),
        dolgo                      = Model('dolgo'),
        _color                     = Model('color'),
        art                        = Model('art'),
        jaeger                     = Model('jaeger'),
        diacritics                 = diacritics,
        model                      = Model('sca'),
        vowels                     = vowels,
        tones                      = tones,
        figsize                    = (10,10),
        combiners                  = '\u0361\u035c',
        breaks                     = '.-',
        stress                     = "ˈˌ'",
        merge_vowels               = True,
        unique_sequences           = True,
        comment                    = '#',
        restricted_chars           = '_T',
        scale                      = 0.5,
        factor                     = 0.3,
        gap_weight                 = 0.5,
        classes                    = True,
        sonar                      = True,
        scorer                     = {},
        tree_calc                  = 'neighbor',
        gop                        = -2,
        ref                        = 'cogid',
        morpheme_separator         = "◦",
        nasal_placeholder          = "∼"
        )
rcParams.update(rcParamsUpd)

# define aliases for parameters
kw_base = dict(
    filename = ('filename', 'fn'),
    merge_vowels = ('mv',),
    sca = ("model",),
    )
alias = {}
for key in kw_base:
    # set key as key, just to make sure that the keyword always occurs in the
    # alias dictionary
    alias[key] = key

    # set all the alias values
    for value in kw_base[key]:
        alias[value] = key

# function changes parameters
def rc(rval=None, **keywords):
    """
    Function changes parameters globally set for LingPy sessions.

    Parameters
    ----------
    rval : string (default=None)
        Use this keyword to specify a return-value for the rc-function. 
    schema : {"ipa", "asjp"}
        Change the basic schema for sequence comparison. When switching to
        "asjp", this means that sequences will be treated as sequences in ASJP
        code, otherwise, they will be treated as sequences written in basic
        IPA.
    verbose : bool (default=False)
        Use this keyword in order to switch to verbose output. This will be
        useful when using complex methods, in order to understand what the
        program is actually doing.
    debug : bool (default=False)
        Use this keyword to switch to debug-mode. It will give specific,
        internal output that is much more technical than the output resulting
        from "verbose".

    Notes
    -----
    This function is the standard way to communicate with the *rcParams*
    dictionary which is not imported as a default. If you want to see which
    parameters there are, you can load the rcParams dictonary directly::

    >>> from lingpy.settings import rcParams

    However, be careful when changing the values. They might produce some
    unexpected behavior.

    Examples
    --------
    Import LingPy:

    >>> from lingpy import *

    Change basic values. Switch to verbose output, for example:

    >>> rc(verbose=True)
    [i] Successfully changed parameters.
    
    """
    from . import log

    if rval:
        return rcParams[rval]
    
    for key in keywords:
        if key == "schema":
            if keywords[key] in ["qlc",'ipa']:
                diacritics,vowels,tones = load_dvt(path='')
                rcParams['asjp'] = Model('asjp')
                rcParams['sca'] = Model('sca')
                rcParams['dolgo'] = Model('dolgo')
                rcParams['art'] = Model('art')
                rcParams['diacritics'] = diacritics
                rcParams['vowels'] = vowels
                rcParams['tones'] = tones
                rcParams['_color'] = Model('color')
                rcParams['combiners']    = '\u0361\u035c'
                rcParams['breaks']       = '.-'
                rcParams['stress']       = "ˈˌ'"
                rcParams['merge_vowels'] = True
                rcParams['basic_orthography'] = 'fuzzy'

                # reset basic model to sca
                rcParams['model'] = rcParams['sca']
                
            elif keywords[key] in ['evolaemp','el','asjp']:
                diacritics,vowels,tones = load_dvt(path='el')
                rcParams['asjp'] = Model('asjp_el')
                rcParams['sca'] = Model('sca_el')
                rcParams['dolgo'] = Model('dolgo_el')
                rcParams['art'] = Model('art_el')
                rcParams['jaeger'] = Model('jaeger_el')
                rcParams['diacritics'] = diacritics
                rcParams['vowels'] = vowels
                rcParams['tones'] = tones
                rcParams['_color'] = Model('color_el')
                rcParams['combiners']    = '\u0361\u035c'
                rcParams['breaks']       = '.-'
                rcParams['stress']       = "ˈˌ'"
                rcParams['merge_vowels'] = False
                rcParams['basic_orthography'] = 'asjp'

                # reset the basic model to the asjp model
                rcParams['model'] = rcParams['asjp']

        if key in alias:
            rcParams[alias[key]] = keywords[key]
        else:
            rcParams[key] = keywords[key]
    log.info("Successfully changed parameters.")
