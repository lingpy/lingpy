"""
Example script for the spreadsheet reader.
"""

__author__ = "Steven Moran"
__date__ = "2014-28-02"

import sys
from lingpy.basic import *

# set debug and/or verbose mode on or off
rcParams["debug"]=True
# rcParams["verbose"]=True

# load the spreadsheet and specify attributes like the blacklist file and the column for concepts
s = Spreadsheet(
    "dogon_wordlists.tsv", 
    meanings="Leipzig-Jakarta", 
#    meanings="Swadesh (AH)", 
#    meanings="English", 
    skip_empty_concepts=False, 
    cellsep="\\\\", 
    blacklist="dogon.bl")

# load as LingPy wordlist, tokenize and create QLC output format
wl = Wordlist(s)
wl.tokenize("Heath2014", column="ipa")
wl.output('qlc',filename='heath2014-LJ-tokenized')

# get some stats from the spreadsheet
# s.stats()
# analyses = s.analyze("words")
# analyses = s.analyze("graphemes")
# s.pprint(analyses[0])
