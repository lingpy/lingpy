# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-17 10:40
# modified : 2013-07-18 12:46
"""
Module handels all global parameters used in a LingPy session.
"""

__author__="Johann-Mattis List"
__date__="2013-07-18"

# builtin imports
from datetime import datetime,date
import os

# internal imports
from ._settings import rcParams
from .data.model import Model,load_dvt

# load diacritics, vowels, tones
diacritics, vowels, tones = load_dvt()

# set the absolute path to lingpy
_abs_path = os.path.split(
        os.path.abspath(
            __file__
            )
        )[0]

# these are general warnings, all prefixed by "warning" in rcParams
warnings = dict(
        warning_empty_cons = "[WARNING] There are emtpy segments in the consensus.",
        warning_failed_cons = "[WARNING] Failed to compute the consensus string.",
        warning_deprecation = "[WARNING] Use of '{0}' is deprecated, use '{1}' instead.",
        warning_missing_module = "[WARNING] Module '{0}' could not be loaded.  Some methods may not work properly.",
        warning_identical_scorer = "[WARNING] An identical scoring function has already been calculated, force recalculation by setting 'force' to 'True'.",
        warning_overwrite_scorer = "[WARNING] A different scoring function has already been calculated, overwriting previous settings.",
        )
rcParams.update(warnings)

# these are questions which need to be answered by the user
questions = dict(
        )
rcParams.update(questions)

# these are lexstat-specific parameters, all prefixed by "lexstat"
lexstat = dict(
        lexstat_transform =  {                    
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
        lexstat_runs = 1000,
        lexstat_modes = [("global",-2,0.5),("local",-1,0.5)],        
        lexstat_rands = 1000,
        lexstat_limit = 10000,
        lexstat_method = 'markov',
        lexstat_ratio = (2,1),
        lexstat_vscale = 0.5,
        lexstat_threshold = 0.7,
        lexstat_cluster_method = 'upgma',
        )
rcParams.update(lexstat)

# these are alignment-specific parameters, all prefixed by "align"
alignments = dict(
        align_mode                       = 'global',
        align_modes                      = [
                ('global',-2,0.5),
                ('local',-1,0.5),
                ],
        align_scale                      = 0.5,
        align_factor                     = 0.3,
        align_gap_weight                 = 0.5,
        align_classes                    = True,
        align_sonar                      = True,
        align_scorer                     = {},
        align_tree_calc                  = 'neighbor',
        )
rcParams.update(alignments)

# dictionary stores basic parameters that are used during a LingPy session
rcParamsUpd = dict(
        verbose                    = False,
        schema                     = 'qlc',
        asjp                       = Model('asjp'),
        sca                        = Model('sca'),
        dolgo                      = Model('dolgo'),
        _color                     = Model('color'),
        art                        = Model('art'),
        diacritics                 = diacritics,
        vowels                     = vowels,
        tones                      = tones,
        figsize                    = (10,10),
        file_written               = "[i] Data has been written to file <{0}>.",
        missing_module             = "[WARNING] Module '{0}' could not be loaded. Some methods may not work properly.",
        filename                   = 'lingpy-'+str(date.today()),
        deprecation_warning        = "[WARNING] Use of '{0}' is deprecated, use '{1}' instead.",
        timestamp                  = str(datetime.today()),
        combiners                  = '\u0361\u035c',
        breaks                     = '.-',
        stress                     = "ˈˌ'",
        merge_vowels               = True,
        _path                      = _abs_path,
        comment                    = '#',
        restricted_chars           = '_T',
        align_mode                       = 'global',
        align_modes                      = [
                ('global',-2,0.5),
                ('local',-1,0.5),
                ],
        scale                      = 0.5,
        factor                     = 0.3,
        gap_weight                 = 0.5,
        classes                    = True,
        sonar                      = True,
        scorer                     = {},
        tree_calc                  = 'neighbor',
        gop                        = -2,
        ref                        = 'cogid',
        #lexstat_transform =  {                    
        #            'A':'C', 
        #            'B':'C',
        #            'C':'C',
        #            'L':'c',
        #            'M':'c',
        #            'N':'c',
        #            'X':'V', #
        #            'Y':'V', #
        #            'Z':'V', #
        #            'T':'T', #
        #            '_':'_'
        #            },
        #lexstat_runs = 1000,
        #lexstat_modes = [("global",-2,0.5),("local",-1,0.5)],        
        #lexstat_rands = 1000,
        #lexstat_limit = 10000,
        #lexstat_method = 'markov',
        #lexstat_ratio = (2,1),
        #lexstat_vscale = 0.5,
        #lexstat_threshold = 0.7,
        #lexstat_cluster_method = 'upgma',
        identical_scorer_warning = "[i] An identical scoring function has already been calculated, force recalculation by setting 'force' to 'True'.",
        overwrite_scoring_function = "[i] A different scoring function has already been calculated, overwriting previous settings.",
        empty_consensus_warning = '[WARNING] There are empty segments in the consensus!',
        sonority_consensus_warning = '[WARNING] Sonority profile consensus could not be calculated!',
        debug = False
        )

rcParams.update(rcParamsUpd)

# define aliases for parameters
kw_base = dict(
    filename = ('filename','fn'),
    file_written = ('filewritten','fw','file_written'),
    merge_vowels = ('mv','merge_vowels'),
    sca = ("sca","model")
    )
alias = {}
for key in kw_base:
    # set key as key, just to make sure that the keyword always occurs in the
    # alias dictionary
    alias[key] = key

    # set all the alias values
    for value in kw_base[key]:
        alias[value] = key

# apply aliases to initial rcParams
for key in list(alias.keys()):
    if key not in rcParams:
        rcParams[key] = rcParams[alias[key]]


# add specific parameters for phybo



# function changes parameters
def rc(**keywords):
    """
    Function changes parameters globally set for LingPy sessions.
    """
    
    for key in keywords:
        if key in rcParams:
            # check for special keyword "schema"
            if key == "schema":
                if keywords[key] == "qlc":
                    diacritis,vowels,tones = load_dvt()
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
                elif keywords[key] in ['evolaemp','el']:
                    diacritics,vowels,tones = load_dvt(path='el')
                    rcParams['asjp'] = Model('asjp_el')
                    rcParams['sca'] = Model('sca_el')
                    rcParams['dolgo'] = Model('dolgo_el')
                    rcParams['art'] = Model('art_el')
                    rcParams['diacritics'] = diacritics
                    rcParams['vowels'] = vowels
                    rcParams['tones'] = tones
                    rcParams['_color'] = Model('color_el')
                    rcParams['combiners']    = '\u0361\u035c'
                    rcParams['breaks']       = '.-'
                    rcParams['stress']       = "ˈˌ'"
                    rcParams['merge_vowels'] = False

            if key in alias:
                for k in alias[key]:
                    rcParams[k] = keywords[key]
            else:
                rcParams[key] = keywords[key]
    if rcParams['verbose']:
        print("[i] Successfully changed parameters.")


