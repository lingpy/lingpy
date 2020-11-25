"""
Module handels all global parameters used in a LingPy session.
"""
from lingpy._settings import rcParams
from lingpy.data.model import Model, load_dvt

# load diacritics, vowels, tones
diacritics, vowels, tones = load_dvt()

# these are lexstat-specific parameters, all prefixed by "lexstat"
lexstat = dict(
    lexstat_transform={
        'A': 'C',
        'B': 'C',
        'C': 'C',
        'L': 'c',
        'M': 'c',
        'N': 'c',
        'X': 'V',
        'Y': 'V',
        'Z': 'V',
        'T': 'T',
        '_': '_'
    },
    lexstat_runs=1000,
    lexstat_modes=[("global", -2, 0.5), ("local", -1, 0.5)],
    lexstat_rands=1000,
    lexstat_limit=10000,
    lexstat_scoring_method='shuffle',
    lexstat_ratio=(2, 1),
    lexstat_vscale=1.0,
    lexstat_threshold=0.45,
    lexstat_cluster_method='upgma',
    lexstat_preprocessing_method='sca',
    lexstat_preprocessing_threshold=0.7,
    lexstat_bad_chars_limit=0.1,
    lexstat_scoring_threshold=0.7
)
rcParams.update(lexstat)

# these are alignment-specific parameters, all prefixed by "align"
alignments = dict(
    align_mode='global',
    align_modes=[('global', -2, 0.5), ('local', -1, 0.5)],
    align_scale=0.5,
    align_factor=0.3,
    align_gap_weight=0.5,
    align_classes=True,
    align_sonar=True,
    align_scorer={},
    align_tree_calc='neighbor',
    align_gop=-2,
    align_transform={
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
    },
    align_notransform={
        # new values for alternative prostrings
        'A': 1,  # initial
        'B': 1,  # syllable-initial
        'C': 1,  # ascending
        'L': 1,  # descending
        'M': 1,  # syllable-descending
        'N': 1,  # final
        'X': 1,  # vowel in initial syllable
        'Y': 1,  # vowel in non-final syllable
        'Z': 1,  # vowel in final syllable
        'T': 1,  # Tone
        '_': 1  # break character
    },
    align_stamp="""# MSA
# dataset    : {0}
# collection : {1}
# aligned by : LingPy Version {2} <www.lingpy.org>
# created on : {3}
# parameters : {4}
""")
rcParams.update(alignments)

# dictionary stores basic parameters that are used during a LingPy session
rcParamsUpd = dict(
    schema='qlc',
    asjp=Model('asjp'),
    sca=Model('sca'),
    dolgo=Model('dolgo'),
    _color=Model('color'),
    art=Model('art'),
    cv=Model('cv'),
    jaeger=Model('jaeger'),
    diacritics=diacritics,
    model=Model('sca'),
    vowels=vowels,
    tones=tones,
    figsize=(10, 10),
    combiners='\u0361\u035c',
    breaks='.-',
    stress="ˈˌ'",
    merge_vowels=True,
    unique_sequences=True,
    comment='#',
    restricted_chars='_T',
    scale=0.5,
    factor=0.3,
    gap_weight=0.5,
    classes=True,
    sonar=True,
    scorer={},
    tree_calc='neighbor',
    gop=-2,
    ref='cogid',
    morpheme_separator="+",
    morpheme_separators="◦+→←",
    nasal_placeholder="∼",
    gap_symbol="-",
    internal_morpheme_separator='_',
    word_separator="_",
    word_separators="_#",
)
rcParams.update(rcParamsUpd)

# define aliases for parameters
kw_base = dict(
    filename=('filename', 'fn'),
    merge_vowels=('mv',),
    sca=("model",),
)
alias = {}
for key in kw_base:
    # set key as key, just to make sure that the keyword always occurs in the
    # alias dictionary
    alias[key] = key

    # set all the alias values
    for value in kw_base[key]:
        alias[value] = key


def rc(rval=None, rcParams_=None, **keywords):
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
    rcParams_ : Allow passing in a plain `dict` for testing.

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

    Switch from IPA transcriptions to ASJP transcriptions:

    >>> rc(schema="asjp")

    You can check which "basic orthography" is currently loaded:

    >>> rc(basic_orthography)
    'asjp'
    >>> rc(schema='ipa')
    >>> rc(basic_orthography)
    'fuzzy'

    """
    from lingpy import log

    # By default we change the global `lingpy._settings.rcParams`
    rcParams_ = rcParams if rcParams_ is None else rcParams_
    if rval:
        return rcParams_[rval]

    for key in keywords:
        if key == "schema":
            if keywords[key] in ["qlc", 'ipa']:
                diacritics, vowels, tones = load_dvt(path='')
                rcParams_['asjp'] = Model('asjp')
                rcParams_['sca'] = Model('sca')
                rcParams_['dolgo'] = Model('dolgo')
                rcParams_['art'] = Model('art')
                rcParams_['diacritics'] = diacritics
                rcParams_['vowels'] = vowels
                rcParams_['tones'] = tones
                rcParams_['_color'] = Model('color')
                rcParams_['combiners'] = '\u0361\u035c'
                rcParams_['breaks'] = '.-'
                rcParams_['stress'] = "ˈˌ'"
                rcParams_['merge_vowels'] = True
                rcParams_['basic_orthography'] = 'fuzzy'

                # reset basic model to sca
                rcParams_['model'] = rcParams['sca']

            elif keywords[key] in ['evolaemp', 'el', 'asjp']:
                diacritics, vowels, tones = load_dvt(path='el')
                rcParams_['asjp'] = Model('asjp_el')
                rcParams_['sca'] = Model('sca_el')
                rcParams_['dolgo'] = Model('dolgo_el')
                rcParams_['art'] = Model('art_el')
                rcParams_['jaeger'] = Model('jaeger_el')
                rcParams_['diacritics'] = diacritics
                rcParams_['vowels'] = vowels
                rcParams_['tones'] = tones
                rcParams_['_color'] = Model('color_el')
                rcParams_['combiners'] = '\u0361\u035c'
                rcParams_['breaks'] = '.-'
                rcParams_['stress'] = "ˈˌ'"
                rcParams_['merge_vowels'] = False
                rcParams_['basic_orthography'] = 'asjp'

                # reset the basic model to the asjp model
                rcParams_['model'] = rcParams['asjp']

        if key in alias:
            rcParams_[alias[key]] = keywords[key]
        else:
            rcParams_[key] = keywords[key]
    log.info("Successfully changed parameters.")
