"""
QLC wordlist module for orthographic parsing.
"""

__author__="Steven Moran"
__date__="2013-03"

# don't use absolute imports!!! XXX
from .wordlist import *
from ..sequence.orthography import *
from ..sequence.sound_classes import *

import os

class QLCWordlist(Wordlist):
    """
    Extended Wordlist class to tokenize QLC data.

    Parameters
    ----------
    filename : { str dict }
        The input file that contains the data. Otherwise a dictionary with
        consecutive integers as keys and lists as values with the key 0
        specifying the header.

    row : str (default = "concept")
        A string indicating the name of the row that shall be taken as the
        basis for the tabular representation of the word list.
    
    col : str (default = "doculect")
        A string indicating the name of the column that shall be taken as the
        basis for the tabular representation of the word list.
    
    conf : string (default='')
        A string defining the path to the configuration file. 
    
    ortho_profile : str (default='')
        An orthographic profile used to convert and tokenize the input data
        into IPA tokens.

    Notes
    -----
    This is a direct daughter of the
    :py:class:`~lingpy.basic.wordlist.Wordlist` class. The difference to a
    "normal" wordlist is that this class searches for an orthography profile
    first, and, if this is defined, tokenizes the input data and adds a new
    column to the wordlist.

    """
    def __init__(
            self, 
            filename, 
            row = 'concept',
            col = 'language',
            conf = '',
            ortho_profile=None, 
            verbose = False, # [JML] added verbose option
            **keywords
            ):

        # add default values
        defaults = {
                "ipa_parse" : "op_tokens",
                "grapheme_parse" : "grapheme_tokens"
                }

        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        # initialize the wordlist
        Wordlist.__init__(self,filename)
        
        # path to orthography profiles
        data_path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + '/data/orthography_profiles/'

        if ortho_profile:
            
            # this is just a very quick fix, that is not very elegant, but it
            # suffices for the current purpose, so it should be cleaned up
            # later on [JML]
            if os.path.exists(ortho_profile):
                filenamepath = ortho_profile
            else:
                filenamepath = data_path+ortho_profile

            if os.path.exists(filenamepath):
                if verbose: print("\n[i] Orthography profile exists. Proceeding with orthographic parsing.\n")
                o = OrthographyParser(filenamepath)

                self.add_entries(keywords['grapheme_parse'], "counterpart", lambda x:o.parse_graphemes(x))

                if o.exists_multiple_columns():
                    self.add_entries(keywords['ipa_parse'], "counterpart", lambda x:o.graphemes_to_ipa(x))                
                    
            else:
                if verbose: print("\n[i] No orthography profile exists in the specified location and file name.")
                if vebose: print("\n[i] Proceeding with graphemic parsing.")

        # if no file at location exists, throw a warning and proceed with Unicode python regex grapheme parsing
        else:
            print("\n[i] Warning no orthography profile was given. Proceeding with Unicode grapheme parsing.")
            g = GraphemeParser()
            self.add_entries("qlc_unicode_grapheme_parse", "counterpart", lambda x:g.parse_graphemes(x))
 
