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

# try to import the precompiled models
# XXX split this up into several parts in order to avoid messing around with
# abundant compilations of stuff that has already been compiled
try:
    ipa_diacritics,ipa_vowels,ipa_tones, = load_dvt()
    sca = Model('sca')
    asjp = Model('asjp')
    dolgo = Model('dolgo')
    art = Model('art')
    _color = Model('color')
# compile the models if they are not precompiled
except:
    from .derive import *
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



