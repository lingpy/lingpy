"""
LingPy comes along with many different kinds of predefined data.  When loading
the library, the following data are automatically loaded and can be used in all
applications:

    .. py:data:: ~lingpy.ipa_diacritics : str
    
       The default string of IPA diacritics which is used for the
       tokenization of IPA strings.
    
    .. py:data:: ipa_vowels : str
    
       The default string of IPA vowels which is used for the tokenization of IPA
       strings.
    
    .. py:data:: sca : Model
       
       The SCA sound-class :py:class:`~lingpy.data.model.Model` (see
       :evobib:`List2012`).
    
    .. py:data:: dolgo : Model
    
       The DOLGO sound-class :py:class:`~lingpy.data.model.Model` (see
       :evobib:`Dolgopolsky1986`).
    
    .. py:data:: asjp : Model
    
       The ASJP sound-class :py:class:`~lingpy.data.model.Model` (see
       :evobib:`Brown2008` and :evobib:`Brown2011`).
    
    .. py:data:: art : Model
    
       The ART sound-class :py:class:`~lingpy.data.model.Model` which is
       used for the calculation of sonority profiles and prosodic strings (see
       :evobib:`List2012`).

"""

from .model import *
try:
    from .derive import *
except ImportError:
    print("[!] Not all modules could be loaded. Some functions might not work.")

# try to import the precompiled models
# XXX split this up into several parts in order to avoid messing around with
# abundant compilations of stuff that has already been compiled
#try:
#    ipa_diacritics,ipa_vowels,ipa_tones, = load_dvt()
#    sca = Model('sca')
#    asjp = Model('asjp')
#    dolgo = Model('dolgo')
#    art = Model('art')
#    _color = Model('color')
## compile the models if they are not precompiled
#except:
#    from .derive import *
#    compile_dvt()
#    compile_model('sca')
#    compile_model('dolgo')
#    compile_model('art')
#    compile_model('color')
#    ipa_diacritics,ipa_vowels,ipa_tones = load_dvt()
#    sca = Model('sca')
#    asjp = Model('asjp')
#    dolgo = Model('dolgo')
#    art = Model('art')
#    _color = Model('color')


def set_global_model(model):
    """
    Define the global sound-class models used for the current LingPy session.

    Parameters
    ----------
    model : string {'qlc','evolamp'}
        Select between 'qlc' as the standard model of the QLC group and
        'evolamp' as the standard model of the EvoLamp group.
        
    """
    global ipa_diacritics
    global ipa_vowels
    global ipa_tones
    global sca
    global asjp
    global dolgo
    global art
    global _color

    if model in ['default','standard','qlc']:
        try:
            ipa_diacritics,ipa_vowels,ipa_tones, = load_dvt()
            sca = Model('sca')
            asjp = Model('asjp')
            dolgo = Model('dolgo')
            art = Model('art')
            _color = Model('color')
        # compile the models if they are not precompiled
        except:
            compile_dvt()
            compile_model('sca')
            compile_model('dolgo')
            compile_model('art')
            compile_model('color')
            ipa_diacritics,ipa_vowels,ipa_tones = load_dvt()
            sca = Model('sca')
            asjp = Model('asjp')
            dolgo = Model('dolgo')
            art = Model('art')
            _color = Model('color')
    elif model in ['evolamp']:
        try:
            ipa_diacritics,ipa_vowels,ipa_tones, = load_dvt()
            sca = Model('sca')
            asjp = Model('asjp_internal')
            dolgo = Model('dolgo_internal')
            art = Model('art')
            _color = Model('color')
        # compile the models if they are not precompiled
        except:
            compile_dvt()
            compile_model('sca')
            compile_model('dolgo')
            compile_model('art')
            compile_model('color')
            ipa_diacritics,ipa_vowels,ipa_tones = load_dvt()
            sca = Model('sca')
            asjp = Model('asjp')
            dolgo = Model('dolgo')
            art = Model('art')
            _color = Model('color')

set_global_model('default')
