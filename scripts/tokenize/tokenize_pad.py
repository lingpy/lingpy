# -*- coding: utf-8 -*-
"""
Test script for the orthography / tokenization module.
"""

__author__ = "Steven Moran"
__date__ = "2010-11-06"

import codecs

from lingpy import *
from lingpy.sequence.tokenizer import *

rc(debug=True)

# make sure the pad data is in this directory
infile = codecs.open("pad_data_qlc.qlc", "r", "utf-8")
header = infile.readline() # skip pad header

t = Tokenizer("pad_orthography_profile")

print()
print("ID"+"\t"+"ORIGINAL"+"\t"+"RULES")
for line in infile:
    line = line.strip()
    tokens = line.split("\t")
    id = tokens[0]
    counterpart = tokens[2]
    grapheme_clusters =t.grapheme_clusters(counterpart)
    rules = t.rules(grapheme_clusters)
    modifiers = t.combine_modifiers(rules)
    # print(id+"\t"+counterpart+"\t"+rules)
    print(id+"\t"+counterpart+"\t"+modifiers)


"""
# this tokenize does not work because of the way the orthography rules are currently written, i.e. 
# they expect space delimited tokens; the wordlist.tokenize() function first apples the rules
# and the the Unicode grapheme cluster tokenization

wl = Wordlist("pad_data_qlc.qlc")
wl.tokenize('pad_orthography_profile')
wl.output("qlc", filename="tokenized-pad_data")
"""
