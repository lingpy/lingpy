"""
QLC wordlist module for orthographic parsing.
"""

__author__="Steven Moran"
__date__="2013-03"

from lingpy.basic.wordlist import *
from lingpy.sequence.orthography import *
from lingpy.sequence.sound_classes import *

import os

class QLCWordlist(Wordlist):
    """
    Extended Wordlist class to tokenize QLC data.

    takes the IPA line and an orthography profile and then basically adds a column to the extended wordlist class

    """
    def __init__(self, filename, ortho_profile=None, **keywords):

        # initialize the wordlist
        Wordlist.__init__(self,filename)

        # path to orthography profiles
        data_path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + '/data/orthography_profiles/'

        if ortho_profile:
            filenamepath = data_path+ortho_profile

            if os.path.exists(filenamepath):
                print("\n[i] Orthography profile exists. Proceeding with orthographic parsing.\n")
                o = OrthographyParser(filenamepath)

                self.add_entries("qlc_grapheme_parse", "counterpart", lambda x:o.parse_graphemes(x))

                if o.exists_multiple_columns():
                    self.add_entries("qlc_ipa_parse", "counterpart", lambda x:o.graphemes_to_ipa(x))                
                    
            else:
                print("\n[i] No orthography profile exists in the specified location and file name.")
                print("\n[i] Proceeding with graphemic parsing.")

        # if no file at location exists, throw a warning and proceed with Unicode python regex grapheme parsing
        else:
            print("\n[i] Warning no orthography profile was given. Proceeding with Unicode grapheme parsing.")
            g = GraphemeParser()
            self.add_entries("qlc_unicode_grapheme_parse", "counterpart", lambda x:g.parse_graphemes(x))


if __name__=="__main__":
    qlc_data_file = QLCWordlist("/Users/stiv/Dropbox/Projects/lingpy/scripts/data/huber1992-test.csv", "huber1992.txt")
    # qlc_data_file = QLCWordlist("/Users/stiv/Dropbox/Projects/lingpy/scripts/data/huber1992-test.csv")
    # print(qlc.header)
    qlc_data_file.output('csv')    


