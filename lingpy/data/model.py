"""
Module for handling sequence models.
"""
from __future__ import division, print_function, unicode_literals
import os
import re

from lingpy.data.derive import compile_model, compile_dvt
from lingpy.read import read_scorer
from lingpy import cache
from lingpy import compat
from lingpy import util


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
    are automatically loaded, and they are accessible via the
    :py:func:`~lingpy.settings.rc` function for
    global settings:

    >>> from lingpy import *
    >>> rc('asjp')
    <sca-model "asjp">

    Define variables for the standard models for convenience:
    
    >>> asjp = rc('asjp')
    >>> sca = rc('sca')
    >>> dolgo = rc('dolgo')
    >>> art = rc('art')

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
    lingpy.data.derive.compile_dvt

    """

    def __init__(self, model, path=None):
        new_path = lambda *cmps: \
            os.path.join(path or util.data_path('models'), model, *cmps)
        self.name = model

        # try to load the converter
        try:
            self.converter = cache.load(model + '.converter')
        except:
            compile_model(model, path)
            self.converter = cache.load(model + '.converter')

        # give always preference to scorer matrix files
        if os.path.isfile(new_path('matrix')):
            self.scorer = read_scorer(new_path('matrix'))
        elif os.path.isfile(new_path('scorer.bin')):
            try:
                self.scorer = cache.load(model + '.scorer')
            except compat.FileNotFoundError:
                pass
        # if none of the above fits, leave it
        else:
            pass

        # read information from the info-file
        self.info = {}

        info = util.read_text_file(new_path('INFO'))
        data = ['description', 'compiler', 'source', 'date', 'vowels', 'tones']

        for line in data:
            try:
                self.info[line] = re.findall('@' + line + ': (.*)', info)[0]
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
        return '<sca-model "' + self.name + '">'

    def __getitem__(self, x):
        return self.converter[x]

    def __contains__(self, x):
        return x in self.converter

    def __eq__(self, x):
        """
        Compare a sound-class model with another model.
        """
        return self.__repr__() == x.__repr__()

    def __call__(self, x, y):
        """
        Use the call-shortcut to retrieve the scoring function.
        """
        return self.scorer[x, y]


def load_dvt(path=''):
    """
    Function loads the default characters for IPA diacritics and IPA vowels of LingPy.
    """
    # check for specific evolaemp path which sues asjp alphabet instead of IPA
    if path in ['el', 'evolaemp']:
        fn = 'dvt_el'
    else:
        fn = 'dvt'

    try:
        dvt = cache.load(fn)
    except:
        compile_dvt(path)
        dvt = cache.load(fn)

    return dvt
