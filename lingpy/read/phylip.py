# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-05 08:22
# modified : 2013-08-19 17:31
"""
Module provides functions to read in various formats from the Phylip package.
"""

__author__="Johann-Mattis List"
__date__="2013-08-19"

from ..settings import *

try:
    import regex as re
except ImportError:
    import re
    print(rcParams['W_missing_module'].format('regex'))
    
import os

from ..algorithm import misc
from .csv import csv2list

def read_dst(
        filename
        ):
    """
    Function reads files in Phylip dst-format.

    Parameters
    ----------
    filename : string
        Name of the file which should have the extension ``dst``.

    Returns
    -------
    data : tuple
        A tuple consisting of a list of taxa and a matrix.

    """
    if os.path.isfile(filename):
        f = open(filename)
    else:
        f = filename.split('\n') # XXX temporary solution
        #print("[!] Could not find the file {0}!".format(filename))

    taxa,matrix = [],[]
    
    
    for i,line in enumerate(f):
        if i > 0:
            taxa.append(line[0:10].strip())
            matrix.append([float(i) for i in
                re.split('\s+',line[11:].strip())])

    return taxa,matrix

def read_scorer(infile):
    """
    Read a scoring function in a file into a ScoreDict object.
    """
    # XXX note that we need a better check here, the previous version also
    # caused a very hard-to-track bug in our system! XXX
    if "\t" in infile and "\n" in infile:
        data = [x.split('\t') for x in infile.split('\n') if x]
    else:
        # read data
        data = csv2list(infile)

    # get the chars
    chars = [l[0] for l in data]
    
    # get the matrix
    matrix = [[float(x) for x in l[1:]] for l in data if l]

    return misc.ScoreDict(chars,matrix)

