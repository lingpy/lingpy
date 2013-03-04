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
words = wl.get_list(col='Russian',flat=True,entry='IPA')

# tokenize all words
for istring in words:
    tstring = ipa2tokens(istring)
    cstring = tokens2class(tstring,sca)
    pstring = prosodic_string(tstring)

    print(
            '{0:10}\t{1:10}\t{2:20}\t{3:10}'.format(
                ''.join(istring),
                ''.join(cstring),
                ' '.join(tstring),
                ''.join(pstring)
                )
            )

for istring in words:
    tstring = ipa2tokens(istring)
    pstring = prosodic_string(tstring)
    sstring = tokens2class(tstring,art)
    if '?' in pstring:
        print('\t'.join(tstring))
        print('\t'.join(pstring))
        print('\t'.join(sstring))
        print('')


tstring = ipa2tokens(argv[1])
pstring = prosodic_string(tstring)
sstring = tokens2class(tstring,art)
print('\t'.join(tstring))
print('\t'.join(pstring))
print('\t'.join(sstring))
print('')

