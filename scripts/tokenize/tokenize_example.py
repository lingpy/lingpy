# -*- coding: utf-8 -*-
"""
Test script for the orthography / tokenization module.
"""

__author__ = "Steven Moran"
__date__ = "2010-12-01"

from lingpy import *

word = "Màttís List"

# t = Tokenizer("../../data/orthography_profiles/test.prf")
# t = Tokenizer("test.prf")
t = Tokenizer("test2")

print()
print(word+"\t\t\t"+"word")
print(t.characters(word)+"\t\t"+"Tokenizer.characters(word)")
print(t.grapheme_clusters(word)+"\t\t"+"Tokenizer.grapheme_clusters(word)")
print(t.graphemes(word)+"\t\t"+"Tokenizer.graphemes(word)")
print(t.transform(word)+"\t\t"+"Tokenizer.transform(word) # defaults to graphemes")
print(t.transform(word, "ipa")+"\t\t"+"Tokenizer.transform(word, 'ipa')")
print(t.transform(word, "funny")+"\t\t"+"Tokenizer.transform(word, 'funny')")
print(t.rules(word)+"\t\t\t\t"+"Tokenizer.rules(word)")
print(t.transform_rules(word)+"\t\t\t"+"Tokenizer.transform_rules(word)")

"""
print()
t2 = Tokenizer()
print()
print(word+"\t\t\t"+"word")
print(t2.characters(word)+"\t\t"+"Tokenizer.characters(word)")
print(t2.grapheme_clusters(word)+"\t\t"+"Tokenizer.grapheme_clusters(word)")
"""
