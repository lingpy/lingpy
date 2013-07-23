# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-06 23:14
# modified : 2013-07-17 14:51
"""
Module for handling sequence models.
"""

__author__="Johann-Mattis List"
__date__="2013-07-17"


import re
from pickle import load
import os
import codecs

from .._settings import rcParams
try:
    from .derive import compile_model,compile_dvt
except ImportError:
    print(rcParams['W_missing_module'].format("networkx"))
try:
    from ..algorithm.cython import misc
except:
    from ..algorithm.cython import _misc as misc
from ..read import *
from ..convert import *

class Model(object):
    """
    Class for the handling of sound-class models.
    
    Parameters
    ----------

    model : { 'sca', 'dolgo', 'asjp', 'art', '_color' }
        A string indicating the name of the model which shall be loaded.
        Select between:
        
        * 'sca' - the SCA sound-class model (see :evobib:`List2012a`),
        * 'dolgo' - the DOLGO sound-class model (see:
          :evobib:`Dolgopolsky1986'),
        * 'asjp' - the ASJP sound-class model (see
          :evobib:`Brown2008` and :evobib:`Brown2011`),
        * 'art - the sound-class model which is used for the calculation of
          sonority profiles and prosodic strings (see :evobib:`List2012`), and
        * '_color" - the sound-class model which is used for the coloring of
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
    are automatically loaded, and they are accessible as global variables via
    their names:

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
            path = ""
            ):
        
        if not path:
            new_path = os.path.join(
                    rcParams['_path'],
                    'data',
                    'models',
                    model
                    )
        else:
            new_path = os.path.join(path,model)
        self.name = model
        # check for converter
        if not os.path.isfile(os.path.join(new_path,'converter.bin')):
            compile_model(model,path)
        
        self.converter = load(open(os.path.join(new_path,'converter.bin'),'rb'))

        # check for scorer
        
        # give always preference to scorer matrix files
        if os.path.isfile(os.path.join(new_path,'matrix')):
            self.scorer = read_scorer(os.path.join(new_path,'matrix'))
        elif os.path.isfile(os.path.join(new_path,'scorer.bin')):
            self.scorer = load(open(os.path.join(new_path,'scorer.bin'),'rb'))

        # if none of the above fits, leave it
        else: 
            pass
        
        # read information from the info-file
        self.info = {}
        
        info = codecs.open(os.path.join(new_path,'INFO'),'r','utf-8').read()
        
        data = ['description','compiler','source','date','vowels','tones']
        
        for line in data:
            try:
                self.info[line] = re.findall('@'+line+': (.*)',info)[0]
            except:
                self.info[line] = 'unknown'

        # check for vowels and tones
        if "vowels" in self.info:
            self.vowels = self.info['vowels']
        if "tones" in self.info:
            self.tones = self.info['tones']

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
        pathx = os.path.join(
                rcParams['_path'],
                'data',
                'models',
                'dvt',
                'dvt.bin'
                )
    elif path in ['el','evolaemp']:
        pathx = os.path.join(
                rcParams['_path'],
                'data',
                'models',
                'dvt_el',
                'dvt.bin'
                )
    else:
        pass
    
    try:
        dvt = load(open(pathx,'rb'))
    except:
        compile_dvt(path)
        dvt = load(open(pathx,'rb'))

    return dvt
