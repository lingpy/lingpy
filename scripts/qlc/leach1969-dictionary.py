"""
Example script that loads and tokenizes QLC data into a LingPy dictionary.

1. Load Leach 1969 data from QLC format into a LingPy dictionary object
2. Tokenize the data given the .prf and .rules files in data/orthography_profiles
3. Write the tokenized data to disk in QLC format.

"""

__author__ = "Steven Moran"
__date__ = "2013-11-01"

from lingpy import *
from lingpy.read import read_qlc

f = read_qlc("../../lingpy/data/qlc/leach1969-67-161.csv")
d = Dictionary(f)
d.tokenize("leach1969")
# d.output("qlc", filename="tokenized-leach1969")



