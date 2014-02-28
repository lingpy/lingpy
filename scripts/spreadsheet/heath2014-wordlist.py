"""
Example script for the spreadsheet reader.
"""

__author__ = "Steven Moran"
__date__ = "2014-28-02"

import sys
from lingpy.basic import *

rcParams["debug"]=True
# rcParams["verbose"]=True

# load the spreadsheet and specify attributes like the blacklist file
s = Spreadsheet(
    "dogon-wordlist.tsv", 
    meanings="Leipzig-Jakarta", 
    skip_empty_concepts=False, 
    cellsep="\\\\", 
    blacklist="dogon.bl")

# load as LingPy wordlist, tokenize and create QLC output format
wl = Wordlist(s)
wl.tokenize("Heath2014", column="ipa")
wl.output('qlc',filename='heath2014-tokenized')


# s.stats()
# analyses = s.analyze("graphemes")
# analyses = s.analyze("graphemes")
# analyses = s.analyze("words")
# s.pprint(analyses[0])

# s.pprint_matrix()
# s.transform(**d)
