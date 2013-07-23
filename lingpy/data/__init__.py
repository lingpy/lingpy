"""
LingPy comes along with many different kinds of predefined data.  When loading
the library, the following data are automatically loaded and can be used in all
applications:

    .. py:data:: ipa_diacritics : str
    
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

    .. py:data:: rcParams : dict

       As an alternative to all global variables, this dictionary
       contains all these variables, and additional ones. This dictionary is used
       for internal coding purposes and stores parameters that are globally set (if
       not defined otherwise by the user), such as
          
          * specific debugging messages (warnings, messages, errors)
          * specific flags (verbose, debug)
          * default values, such as "gop" (gap opening penalty), "scale" (scaling factor
          * by which extended gaps are penalized), or "figsize" (the default size of
          * figures if data is plotted using matplotlib).
          
       These default values can be changed with help of the ``rc`` function that takes any
       keyword and any variable as input and adds or modifies the specific key of the
       rcParams dictionary, but also provides more complex functions that change whole
       sets of variables, such as the following statement::
          
          >>> rc(schema="evolaemp")
          
       which switches the variables "asjp", "dolgo", etc. to the ASCII-based
       transcription system of the ASJP project.

"""

from .model import *
