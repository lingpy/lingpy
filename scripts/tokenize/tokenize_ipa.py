"""
Test script for the orthography / tokenization module.
"""

__author__ = "Steven Moran"
__date__ = "2013-11-02"

import codecs

from lingpy.sequence.tokenizer import *

t = Tokenizer()
t.tokenize_ipa("string")

with codecs.open("eng_ladefoged1999-narrow.txt", "r", "utf-8") as file:
    for line in file:
        line = line.strip()
        tokens = line.split()
        for token in tokens:
            token = token.strip()
            print(token+"\t"+t.tokenize_ipa(token))


