# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 13:58
# modified : 2013-03-04 13:58
"""
This script tests various SCA routines.
"""

__author__="Johann-Mattis List"
__date__="2013-03-04"

from lingpy import *
from lingpy.align.pairwise import *
from sys import argv

# load kesslers wordlist
wl = Wordlist('data/IEL.csv')

# extract German words
wordsA = wl.get_dict(col='German',entry='IPA')
wordsB = wl.get_dict(col='Swedish',entry='IPA')

# tokenize all words
#for istring in words:
#    tstring = ipa2tokens(istring)
#    cstring = tokens2class(tstring,sca)
#    pstring = prosodic_string(tstring,'t')
#    #syls = [str(p) for p in prosodic_string(tstring)[1]]
#    prostring = prosodic_string(tstring,'p')
#    proweights = prosodic_weights(prostring)
#
#    print(
#            '{0:10}\t{1:10}\t{2:20}\t{3:10}\t{4:10}'.format(
#                ''.join(istring),
#                ''.join(cstring),
#                ' '.join(tstring),
#                ''.join(pstring),
#                ''.join(prostring)
#                )
#            )
#    print(proweights)

# align words provisionally
seqs = []
for c in wl.concept: 
    if c in wordsA and c in wordsB:
        wA = wordsA[c][0]
        wB = wordsB[c][0]

        seqs += [[wA,wB]]

        # get the tokens
        #tA = ipa2tokens(wA)
        #tB = ipa2tokens(wB)

        #almA,almB,sim = basic_align(tA,tB,distance=True)
        #print(c)
        #print('\t'.join(almA))
        #print('\t'.join(almB))
        #print(sim)
        #print('---')

        ## get the classes
        #cA = tokens2class(tA,sca)
        #cB = tokens2class(tB,sca)

        ## get the prostrings
        #pA = prosodic_string(tA,'p')
        #pB = prosodic_string(tB,'p')

        ## get the weights
        #wgA,wgB = prosodic_weights(pA),prosodic_weights(pB)

        ## align the stuff
        #almA,almB,sim = sc_align(
        #        cA,
        #        cB,
        #        wgA,
        #        wgB,
        #        pA,
        #        pB,
        #        -1,
        #        0.3,
        #        0.5,
        #        sca.scorer,
        #        'T_',
        #        'global',
        #        True
        #        )

        ## convert alignments back to original form
        #outA = class2tokens(tA,almA,local=False)
        #outB = class2tokens(tB,almB,local=False)
        #
        #print(c)
        #print('\t'.join(outA))
        #print('\t'.join(outB))
        #print(sim)
        #print('-----')

from lingpy.align.pairwise import Pairwise as pw
from lingpy.align.multiple import Multiple
#pairs = pw(seqs)
#pairs.align(distance=True,pprint=True)

msa = Multiple(
        [
            'muter',
            'moθər',
            'mur'
            ]
        )

msa._set_model()
msa._set_scorer('classes')
msa._get_pairwise_alignments()
msa._create_library()
msa._extend_library()
msa._make_guide_tree()
msa._merge_alignments()
