# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-09-21 09:06
# modified : 2014-09-21 09:06
"""
Basic parser for Starling data.
"""

__author__="Johann-Mattis List"
__date__="2014-09-21"

import codecs
import os
from ..settings import rcParams
from ..sequence.sound_classes import ipa2tokens
from .csv import csv2list

def star2qlc(filename):
    """
    Converts a file directly output from starling to LingPy-QLC format.
    """

    data = csv2list(filename)

    # get the header
    header = data[0]

    # search for '#' char in header
    cognates = False
    for h in header:
        if '#' in h:
            cognates = True
    # determine language names in header   
    taxa = []
    for i in range(len(header)-1):

        prev = header[i]
        post = header[i+1]
        
        if prev in post and '#' in post:
            taxa += [prev]

            if len(taxa) == 1:
                lngIdx = i

        if prev == 'Number':
            numIdx = i

        if prev == 'Word':
            wrdIdx = i
    
    if rcParams['verbose']:
        print('starling, indices',lngIdx,numIdx,wrdIdx)
        print('starling, taxa:',taxa)

    # start filling in the dictionary
    D = {}
    
    idx = 1
    cognate_counter = 0
    for line in data[2:]:
        
        gloss = line[wrdIdx]
        gnum = line[numIdx]
        
        cognate_sets = []
        for i in range(lngIdx,len(header),2):
            word = line[i]
            
            if '{' in word:
                ipa = word[:word.index('{')].strip()
                ortho = word[word.index('{')+1:word.index('}')].strip()
            else:
                ipa = word
                ortho = word
            
            cogid = int(line[i+1]) 

            if cogid != 0 and word:
            
                cogid = cogid + cognate_counter


                # append cognate sets, essential for raising the counter
                cognate_sets += [int(cogid)]
                
                taxon = header[i]

                D[idx] = [taxon,gloss,gnum,word,ortho,ipa,cogid]
                idx += 1
            

        max_cog = max(cognate_sets)
        cognate_counter = max_cog + 1

    D[0] = ['DOCULECT','CONCEPT','GLOSSID','WORDINSOURCE','ORTHOGRAPHY','IPA','COGID']

    return D


