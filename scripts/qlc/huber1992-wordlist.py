"""
Example script that loads and tokenizes QLC data into a LingPy dictionary.

1. Load Huber 1992 comparative wordlist data from QLC format into a LingPy wordlist
2. Tokenize the data given the .prf and .rules files in data/orthography_profiles
3. Write the tokenized data to disk in QLC format.

"""

__author__ = "Steven Moran"
__date__ = "2013-11-01"

from lingpy import *

file = read_qlc("../../lingpy/data/qlc/huber1992.csv")
wl = Wordlist(file)
wl.tokenize("../../lingpy/data/orthography_profiles/huber1992.prf", "../../lingpy/data/orthography_profiles/huber1992.rules")
wl.output("qlc", filename="tokenized-huber1992")
