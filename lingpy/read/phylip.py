"""
Module provides functions to read in various formats from the Phylip package.
"""
from __future__ import unicode_literals, print_function, division
import re

from lingpy.algorithm import misc
from lingpy.read.csv import csv2list
from lingpy.util import read_text_file


def read_dst(filename, taxlen=10, comment='#'):
    """
    Function reads files in Phylip dst-format.

    Parameters
    ----------
    filename : string
        Name of the file which should have the extension ``dst``.
    taxlen : int (default=10)
        Indicate how long the taxon names are allowed to be in the file from
        which you want to read. The Phylip package
        only allows taxon names consisting of maximally 10 characters (this is
        the default). Other
        packages, however, allow more. If Phylip compatibility is not important
        for you and you just want to allow for as long taxon names as possible,
        set this value to 0 and make sure to use tabstops as separators between
        values in your matrix file.
    comment : str (default = '#')
        The comment character to be used if your file contains additional
        information which should be ignored.

    Returns
    -------
    data : tuple
        A tuple consisting of a list of taxa and a matrix.

    """
    if '\n' in filename:
        lines = [f for f in filename.split('\n') if f.strip()]
    else:
        lines = read_text_file(filename, normalize="NFC", lines=True)

    taxa, matrix = [], []

    for line in lines[1:]:
        if not line.startswith(comment):
            if taxlen > 0:
                taxa.append(line[:taxlen].strip())
                matrix.append([float(val) for val in
                               re.split('\s+', line[taxlen + 1:].strip())])
            else:
                splits = line.split('\t')
                taxa.append(splits[0])
                matrix.append([float(val.strip()) for val in splits[1:]])

    return taxa, matrix


def read_scorer(infile):
    """
    Read a scoring function in a file into a ScoreDict object.

    Parameters
    ----------
    infile : str
        The path to the input file that shall be read as a scoring dictionary.
        The matrix format is a simple csv-file in which the scoring matrix is
        displayed, with negative values indicating high differences between
        sound segments (or sound classes) and positive values indicating high
        similarity. The matrix should be symmetric, columns should be separated
        by tabstops, and the first column should provide the alphabet for which
        the scoring function is defined.

    Returns
    -------
    scoredict : ~lingpy.algorithm.misc.ScoreDict
        A ScoreDict instance which can be directly passed to LingPy's alignment
        functions.
    """
    # XXX note that we need a better check here, the previous version also
    # caused a very hard-to-track bug in our system! XXX
    if "\t" in infile and "\n" in infile:
        data = [x.split('\t') for x in infile.split('\n') if x]
    else:
        data = csv2list(infile)

    return misc.ScoreDict(
        [l[0] for l in data], [[float(x) for x in l[1:]] for l in data if l])
