# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-11 20:01
# modified : 2013-03-11 20:01
"""
This module provides miscellaneous functions which are mostly used internally.
"""

__author__="Johann-Mattis List"
__date__="2013-03-11"


import sys
from re import sub,findall
from numpy import array,sqrt,zeros
from ..data import *
from ..data.ipa.sampa import reXS,xs

def squareform(x):
    """
    A simplified version of the :py:func:`scipy.spatial.distance.squareform` \
    function.

    Parameters
    ----------

    x : :py:class:`numpy.array` or list
        The one-dimensional flat representation of a symmetrix distance matrix.

    Returns
    -------
    matrix : :py:class:`numpy.array`
        The two-dimensional redundant representation of a symmetric distance matrix.

    """

    l = len(x)

    # calculate the length of the square
    s = int(sqrt(2 * l) + 1)
    
    out = zeros((s,s))
    
    k = 0
    for i in range(s):
        for j in range(s):
            if i < j:
                out[i][j] = x[k]
                out[j][i] = x[k]
                k += 1
    return out

def loadtxt(infile):
    """
    Function imitates the :py:func:`numpy.loadtxt` function.

    Parameters
    ----------
    infile : file
        The input file from which the data is read.

    Returns
    -------
    data : list
        A list object which renders the dimensions of the input file.

    """

    f = open(infile)

    data = []
    for line in f:
        if not line.startswith('#'):
            data.append([x.strip() for x in line.strip().split('\t')])

    # check the data for consistency
    start = len(data[0])
    problems = False
    for i,x in enumerate(data):
        if len(x) != start:
            problems = True
            break

    # if there are inconsistent lines, an integrity check of the data has to be
    # carried out, ask the user whether this is wanted, or not
    if problems:
        answer = raw_input('[!] Input file contains errors. Do you want the errors to be listed separatedly (y/n)? ')
        if answer == 'y':
            for i,x in enumerate(data):
                if len(x) != start:
                    print("... {0} columns in line {1}, expected {2} ...\
                            ...".format(
                                len(x),
                                i+1,
                                start
                                )
                            )
        else:
            pass
    
    return data

class LingpyArray(object):
    """
    An extension of the numpy array object which allows the storage of lists in
    two-dimensional arrays.

    Parameters
    ----------

    input_list : list
        The list which shall be converted in an array-like object.

    """

    def __init__(self,input_list):

        
        h = len(input_list)
        w = len(input_list[0])

        self.array = zeros((h,w),dtype='int')

        self.dictionary = {}

        count = 0
        for i,line in enumerate(input_list):
            for j,itm in enumerate(line):
                self.dictionary[count] = input_list[i][j]
                self.array[i][j] = count
                count += 1

    def __getitem__(self,x):

        ind = self.array[x]

        try:
            return self.dictionary[self.array[x]]
        except:
            try:
                return [self.dictionary[i] for i in ind]
            except:
                return [[self.dictionary[i] for i in j] for j in ind]

def check_tokens(tokens):
    """
    Function checks whether tokens are given in a consistent input format.
    """
    
    errors = []

    for i,token in enumerate(tokens):
        
        # check for conversion within the articulation-model
        try:
            art.converter[token]
        except:
            try:
                art.converter[token[0]]
            except:
                errors.append((i,token))
    
    return errors

def transpose(matrix):

    lA = len(matrix)
    lB = len(matrix[0])

    out = [[matrix[i][j] for i in range(lA)] for j in range(lB)]

    return out

