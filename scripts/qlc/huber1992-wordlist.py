"""
Example script that loads and tokenizes QLC data into a LingPy dictionary.

1. Load Huber 1992 comparative wordlist data from QLC format into a LingPy wordlist
2. Tokenize the data given the .prf and .rules files in data/orthography_profiles
3. Write the tokenized data to disk in QLC format.

"""

__author__ = "Steven Moran"
__date__ = "2013-11-01"

from lingpy import *
from lingpy.read import read_qlc

f = read_qlc("../../lingpy/data/qlc/huber1992.csv")
wl = Wordlist(f)
wl.tokenize("huber1992")
wl.output("qlc", filename="tokenized-huber1992")
