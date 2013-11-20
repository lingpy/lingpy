"""
Orthography profile validator. Checks to make sure that all columns 
are of equal or less length than the header.
"""

import os
import codecs

dir = os.listdir(".")

for filename in dir:
    if filename.endswith(".prf"):
        infile = codecs.open(filename, "r")

        print("processing: ", filename)
        header = []
        for line in infile:
            line = line.strip()
            # skip metadata
            if line.startswith("#"):
                continue
            # grab the header
            if line.lower().startswith("graphemes"):
                header = line.split("\t")
                continue
            # process the other lines
            tokens = line.split("\t")
            if len(tokens) > len(header):
                print(len(header), header)
                print(len(tokens), tokens)

