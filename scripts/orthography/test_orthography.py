# -*- coding: utf-8 -*-

"""
Test script for the orthography module.
"""

__author__ = "Steven Moran"
__date__ = "2010-12-01"

from lingpy import *

# replace these with
# t = Tokenizer()
# t.characters()
# t.grapheme_clusters()
# t.slashx() -> calls t.grapheme_clusters)
# t.graphemes()
# t.transform(column)
# t.transform(x,y)
# t.rules()

word = "Màttís List"
fail = "fail"

t = Tokenizer("../../lingpy/data/orthography_profiles/test.prf", "../../lingpy/data/orthography_profiles/test.rules")


print()
print(word+"\t\t\t"+"word")
print(t.characters(word)+"\t\t"+"Tokenizer.characters(word)")
print(t.grapheme_clusters(word)+"\t\t"+"Tokenizer.grapheme_clusters(word)")
print(t.graphemes(word)+"\t\t"+"Tokenizer.graphemes(word)")
print(t.transform(word)+"\t\t"+"Tokenizer.transform(word)")
print(t.transform(word, "ipa")+"\t\t"+"Tokenizer.transform(word, 'ipa')")
print(t.transform(word, "funny")+"\t\t"+"Tokenizer.transform(word, 'funny')")
print(t.rules(word)+"\t\t\t\t"+"Tokenizer.rules(word)")

print()
sys.exit(1)


print()
print(word+"\t\t\t"+"original word")
print(g.tokenize_characters(word)+"\t"+"GraphemeParser.tokenize_characters()")
print(g.tokenize_graphemes(word)+"\t\t"+"GraphemeParser.tokenize_graphemes()")
print()
print(o.tokenize_graphemes(word)+"\t\t"+"OrthographyParser.tokenize_graphemes()")
# print(o.transform_graphemes(word)+"\t\t"+"OrthographyParser.tokenize_graphemes()")


print()
print(fail+"\t\t\t"+"word fails, characters missing in test.prf")
print(o.tokenize_graphemes(fail)+"\t\t"+"OrthographyParser.tokenize_graphemes()")

print()

"""
test_words = ["aa", "aabch", "aa abch", "aa-abch"]
print()

for word in test_words:
    print("--------------------------")
    print("original word:", word)
    print()
    print(g.tokenize_graphemes(word)+"\t\t"+"GraphemeParser.tokenize_graphemes()")
    print(g.tokenize_characters(word)+"\t\t"+"GraphemeParser.tokenize_characters()")
    print()

    print(o.tokenize_graphemes(word)+"\t\t"+"OrthographyParser.tokenize_graphemes()")

    print(o.transform(word, "funny"))


#    print(o.transform_graphemes("a", "funny")) # +"\t\t\t"+"OrthographyParser.transform_graphemes()")
#    print(o.transform_graphemes("-", "IPA")) # +"\t\t\t"+"OrthographyParser.transform_graphemes()")

    print()

    print(rules.tokenize_string(word)+"\t\t\t"+"OrthographyRulesParser.tokenize_string()")
    print()
"""
