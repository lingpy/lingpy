# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-06 23:14
# modified : 2013-03-06 23:14
"""
Module for handling sequence models.
"""

__author__="Johann-Mattis List"
__date__="2013-03-06"


from re import findall
from pickle import load
import os
import codecs

from ..check.exceptions import *

try:
    from .derive import compile_model
except ImportError:
    ThirdPartyModuleError('networkx').warning()    

class Model(object):
    """
    Class for the handling of sound-class models.
    
    Parameters
    ----------

    model : { 'sca', 'dolgo', 'asjp', 'art', 'color' }
        A string indicating the name of the model which shall be loaded.
        Select between:
        
        * 'sca' - the SCA sound-class model (see :evobib:`List2012a`),
        * 'dolgo' - the DOLGO sound-class model (see:
          :evobib:`Dolgopolsky1986'),
        * 'asjp' - the ASJP sound-class model (see
          :evobib:`Brown2008` and :evobib:`Brown2011`),
        * 'art - the sound-class model which is used for the calculation of
          sonority profiles and prosodic strings (see :evobib:`List2012`), and
        * 'color" - the sound-class model which is used for the coloring of
          sound-tokens when creating html-output.  
    
    Notes
    -----
    Models are loaded from binary files which can be found in the
    :file:`data/models/` folder of the LingPy package. A model has two
    essential attributes: 
    
    * :py:attr:`converter` -- a dictionary with IPA-tokens as keys and
      sound-class characters as values, and
    * :py:attr:`scorer` -- a scoring dictionary with tuples of sound-class
      characters as keys and scores (integers or floats) as values.
    
    
    Attributes
    ----------
    converter : dict
        A dictionary with IPA tokens as keys and sound-class characters as
        values.

    scorer : dict
        A scoring dictionary with tuples of sound-class characters as keys and
        similarity scores as values.

    info : dict
        A dictionary storing the key-value pairs defined in the ``INFO``.

    name : str
        The name of the model which is identical with the name of the folder
        from wich the model is loaded.
    
    Examples
    --------
    When loading LingPy, the models ``sca``, ``asjp``, ``dolgo``, and ``art``
    are automatically loaded:

    >>> from lingpy import *

    Check, how the letter ``a`` is converted in the various models:

    >>> for m in [asjp,sca,dolgo,art]: 
    >>> for m in [asjp,sca,dolgo,art]:
    ...     print('{0} > {1} ({2})'.format('a',m.converter['a'],m.name))
    ... 
    a > a (asjp)
    a > A (sca)
    a > V (dolgo)
    a > 7 (art)
    
    Retrieve basic information of a given model:

    >>> print(sca)
    Model:    sca
    Info:     Extended sound class model based on Dolgopolsky (1986)
    Source:   List (2012)
    Compiler: Johann-Mattis List
    Date:     2012-03
    
    See also
    --------
    lingpy.data.derive.compile_model
    lingpy.data.derive.compile_diacritics_and_vowels

    """

    def __init__(
            self,
            model,
            path = None
            ):
        
        if path == None:
            new_path = os.path.split(os.path.abspath(__file__))[0]+'/models/'+model+'/'
        else:
            if path.endswith('/'):
                new_path = path+model+'/'
            else:
                new_path = path + '/' + model + '/'

        self.name = model
        try:
            self.converter = load(open(new_path+'converter.bin','rb'))
            try:
                self.scorer = load(open(new_path+'scorer.bin','rb'))
            except:
                pass
        except:
            compile_model(model,path)
            self.converter = load(open(new_path+'converter.bin','rb'))
            try:
                self.scorer = load(open(new_path+'scorer.bin','rb'))
            except:
                pass

        
        # read information from the info-file
        self.info = {}
        
        info = codecs.open(new_path+'INFO','r','utf-8').read()
        
        data = ['description','compiler','source','date']
        
        for line in data:
            try:
                self.info[line] = findall('@'+line+': (.*)',info)[0]
            except:
                self.info[line] = 'unknown'

    def __str__(self):
        
        out = 'Model:    {0}\nInfo:     {1}\nSource:   {2}\n'
        out += 'Compiler: {3}\nDate:     {4}'
        out = out.format(
                self.name,
                self.info['description'],
                self.info['source'],
                self.info['compiler'],
                self.info['date']
                )
        return out

    def __repr__(self):

        return '<sca-model "'+self.name+'">'

    def __getitem__(self, x):

        return self.converter[x]

    def __contains__(self,x):

        if x in self.converter:
            return True
        else:
            return False

    def __eq__(self,x):
        """
        Compare a sound-class model with another model.
        """

        if self.__repr__() == x.__repr__():
            return True
        else:
            return False

    def __call__(self,x,y):
        """
        Use the call-shortcut to retrieve the scoring function.
        """
        return self.scorer[x,y]

def load_dvt(path=''):
    """
    Function loads the default characters for IPA diacritics and IPA vowels of LingPy.
    """
    if not path:
        path = os.path.split(os.path.abspath(__file__))[0]+'/models/dvt/dvt.bin'
    elif path in ['el','evolaemp']:
        path = os.path.split(os.path.abspath(__file__))[0]+'/models/dvt_el/dvt.bin'
    else:
        pass

    dvt = load(open(path,'rb'))

    return dvt

def set_global_model2(model='qlc'):
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
    elif model in ['evolamp', 'el']:
        try:
            ipa_diacritics,ipa_vowels,ipa_tones, = load_dvt(path='el')
            sca = Model('sca_el')
            asjp = Model('asjp_el')
            dolgo = Model('dolgo_el')
            art = Model('art_el')
            _color = Model('color_el')
        # compile the models if they are not precompiled
        except:
            compile_dvt(path='el')
            compile_model('sca_el')
            compile_model('dolgo_el')
            compile_model('art_el')
            compile_model('color_el')
            ipa_diacritics,ipa_vowels,ipa_tones = load_dvt(path='el')
            sca = Model('sca_el')
            asjp = Model('asjp_el')
            dolgo = Model('dolgo_el')
            art = Model('art_el')
            _color = Model('color_el')


