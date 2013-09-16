# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-15 21:41
# modified : 2013-09-15 21:41
"""
Module provides basic operations on Wordlist-Objects.
"""

__author__="Johann-Mattis List"
__date__="2013-09-15"

def renumber(wordlist,source,target=''):
    """
    Create numerical identifiers from string identifiers.
    """

    # iterate over wordlist and get all source ids
    sources = []
    for k in wordlist:
        sources += [wordlist[k,source]]

    sources = sorted(set(sources))

    # convert to numbers
    targets = list(range(1,len(sources)+1))

    # add to wordlist
    if not target:
        target = source+'id'

    # make converter
    converter = dict(zip(sources,targets))
    
    wordlist.add_entries(target,source,lambda x:converter[x])

    # add stuff to meta
    wordlist._meta[source+'2'+target] = converter

