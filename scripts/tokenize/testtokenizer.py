# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-04 22:23
# modified : 2013-11-04 22:23
"""
Test the tokenization function.
"""

__author__="Johann-Mattis List"
__date__="2013-11-04"

from lingpy import *

words = [
    "tʰɔxtər",
    "faːtər",
    "p͡fyt͡sə"
    ]

for word in words:
    
    tokA = tokenize(word)
    tokB = tokenize(word, orthography='plain')
    clsA = tokenize(word, model='sca')
    clsB = tokenize(word, model='dolgo')
    clsC = tokenize(word, model='asjp')

    print("{0:10}\t{1}".format(
        "fuzzy",
        "\t".join(tokA)
        ))
    print("{0:10}\t{1}".format(
        "plain",
        "\t".join(tokB)
        ))
    print("{0:10}\t{1}".format(
        "SCA",
        "\t".join(clsA)
        ))
    print("{0:10}\t{1}".format(
        "DOLGO",
        "\t".join(clsB)
        ))
    print("{0:10}\t{1}".format(
        "ASJP",
        "\t".join(clsC)
        ))

rc(schema='asjp')
word = "mama*n"
tokA = tokenize(word)
tokB = tokenize(word, orthography='asjp')
clsA = tokenize(word, model='sca')
clsB = tokenize(word, model='dolgo')
clsC = tokenize(word, model='asjp')

print("{0:10}\t{1}".format(
    "fuzzy",
    "\t".join(tokA)
    ))
print("{0:10}\t{1}".format(
    "plain",
    "\t".join(tokB)
    ))
print("{0:10}\t{1}".format(
    "SCA",
    "\t".join(clsA)
    ))
print("{0:10}\t{1}".format(
    "DOLGO",
    "\t".join(clsB)
    ))
print("{0:10}\t{1}".format(
    "ASJP",
    "\t".join(clsC)
    ))

