"""
Simple script to take text as input and out it as Unicode NFC.
"""

import sys
import codecs
import unicodedata

filename = sys.argv[1]

infile = codecs.open("../"+filename, "r", "utf-8")
outfile = open(filename, "w")

for line in infile:
    line = line.strip()
    line = unicodedata.normalize("NFC", line)
    outfile.write(line+"\n")

