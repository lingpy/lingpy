# author   : Steven Moran                                                                                                         
# email    : steve.moran@lmu.de                                                                                                   
# created  : 2013-01-10                                                                                                           

"""
Test script for the orthography module.
"""

__author__ = "Steven Moran"
__date__ = "2010-12-01"

from lingpy import *

o = OrthographyParser("../../lingpy/data/orthography_profiles/test.prf")
g = GraphemeParser()
test_words = ["aa", "aabuu", "uuabaa auubaa"]
print()
for word in test_words:
    print("original word:", word)
    print("parse_string_to_graphemes_string:", o.parse_string_to_graphemes_string(word))
    print("parse_string_to_ipa_string:", o.parse_string_to_ipa_string(word))
    print("parse_string_to_graphemes:", o.parse_string_to_graphemes(word))
    # print("parse_graphemes:", g.parse_graphemes(word))                                                                      
    print()
