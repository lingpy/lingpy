# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-05 08:22
# modified : 2013-03-05 08:22
"""
Module provides functions to read in various formats from the Phylip package.
"""

__author__="Johann-Mattis List"
__date__="2013-03-05"

import regex as re
import os
try:
    from algorithm.cython import misc
except ImportError:
    from ..algorithm.cython import _misc as misc
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
    if os.path.isfile(infile):
        # read data
        data = csv2list(infile)
    else:
        data = [x.split('\t') for x in infile.split('\n') if x]


    # get the chars
    chars = [l[0] for l in data]
    
    # get the matrix
    matrix = [[float(x) for x in l[1:]] for l in data if l]

    return misc.ScoreDict(chars,matrix)

