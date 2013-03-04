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
from sys import argv

# load kesslers wordlist
wl = Wordlist('data/IEL.csv')

# extract German words
words = wl.get_list(col='Urdu',flat=True,entry='IPA')

# tokenize all words
for istring in words:
    tstring = ipa2tokens(istring)
    cstring = tokens2class(tstring,sca)
    pstring = [p[0] for p in prosodic_string(tstring)]
    syls = [str(p[1]) for p in prosodic_string(tstring)]

    print(
            '{0:10}\t{1:10}\t{2:20}\t{3:10}\t{4:10}'.format(
                ''.join(istring),
                ''.join(cstring),
                ' '.join(tstring),
                ''.join(pstring),
                ' '.join(syls)
                )
            )

from lingpy.algorithm.misc import prosodic_string as pstr
i = 0
for istring in words:
    tstring = ipa2tokens(istring)
    pstring = prosodic_string(tstring,output='pstring')
    sstring = tokens2class(tstring,art)
    pstring2 = pstr([int(s) for s in sstring])
    if list(pstring2) != pstring:
        print(i+1)
        i += 1
        print('\t'.join(tstring))
        print('\t'.join(pstring))
        #print('\t'.join(sstring))
        print('--')
        print('\t'.join(pstring2))
        print('')


#tstring = ipa2tokens(argv[1])
#pstring = prosodic_string(tstring)
#syls = [str(s[1]) for s in pstring]
#pstring = [p[0] for p in pstring]
#sstring = tokens2class(tstring,art)
#print('\t'.join(tstring))
#print('\t'.join(pstring))
#print('\t'.join(sstring))
#print('\t'.join(syls))
#print('\t'.join(prosodic_string(tstring,output='pstring')))
#print('')

