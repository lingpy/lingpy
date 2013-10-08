# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-15 21:41
# modified : 2013-09-15 21:41
"""
Module provides basic operations on Wordlist-Objects.
"""

__author__="Johann-Mattis List"
__date__="2013-09-15"

from ..settings import rcParams
from ..convert import *

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

    if rcParams['verbose']: print("[i] Successfully renumbered {0}.".format(source))


def calculate(
        wordlist,
        data,
        taxa = 'taxa',
        concepts = 'concepts',
        ref = 'cogid',
        **keywords
        ):
    """
    Manipulate a wordlist object by adding different kinds of data.
    
    Parameters
    ----------
    data : str
        The type of data that shall be calculated. Currently supports

        * "tree": calculate a reference tree based on shared cognates
        * "dst": get distances between taxa based on shared cognates
        * "cluster": cluster the taxa into groups using different methods


    """
    defaults = dict(
            distances = False,
            tree_calc = "upgma",
            cluster = "upgma",
            force = False,
            threshold = 0.5,
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # get taxa for current calculation
    these_taxa = eval('wordlist.'+taxa)
    
    # calculate distances
    if data in ['distances','dst']:
        wordlist._meta['distances'] = wl2dst(
                wordlist,
                taxa,
                concepts,
                ref,
                **keywords
                )
    elif data in ['tre','tree','nwk']:
        if 'distances' not in wordlist._meta:
            distances = wl2dst(wordlist,taxa,concepts,ref,**keywords)
        else:
            distances = wordlist._meta['distances']
        if 'tree' in wordlist._meta and not keywords['force']:
            print(rcParams['W_force'].format('Reference tree'))
            return
        wordlist._meta['tree'] = matrix2tree(
                distances,
                these_taxa,
                keywords['tree_calc'],
                keywords['distances']
                )

    elif data in ['groups','cluster']:
        if 'distances' not in wordlist._meta:
            distances = wl2dst(wordlist,taxa,concepts,ref,**keywords)
        else:
            distances = wordlist._meta['distances']
        if 'groups' in wordlist._meta and not keywords['force']:
            print(rcParams['W_force'].format('Distance matrix'))
            return
        wordlist._meta['groups'] = matrix2groups(
                keywords['threshold'],
                distances,
                these_taxa
                )

    if rcParams['verbose']: print("[i] Successfully calculated {0}.".format(data))
        

                
