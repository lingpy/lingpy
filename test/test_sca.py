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
    pstring = [p for p in prosodic_string(tstring)[0]]
    syls = [str(p) for p in prosodic_string(tstring)[1]]
    prostring = prosodic_string(tstring,'p')
    proweights = prosodic_weights(prostring)

    print(
            '{0:10}\t{1:10}\t{2:20}\t{3:10}\t{4:10}'.format(
                ''.join(istring),
                ''.join(cstring),
                ' '.join(tstring),
                ''.join(pstring),
                ' '.join(syls)
                )
            )
    print(proweights)

