# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 17:02
# modified : 2013-06-26 18:08
"""
Module provides functions for reading csv-files.
"""

__author__="Johann-Mattis List"
__date__="2013-06-26"

import codecs
import os

def csv2list(
        filename,
        fileformat = '',
        dtype = [],
        comment = '#',
        sep = '\t',
        strip_lines = True
        ):
    """
    Very simple function to get quick access to CSV-files.

    Parameters
    ----------
    filename : str
        Name of the input file.
    fileformat : {None str}
        If not specified the file <filename> will be loaded. Otherwise, the
        fileformat is interpreted as the specific extension of the input file.
    dtype : {list}
        If not specified, all data will be loaded as strings. Otherwise, a
        list specifying the data for each line should be provided.
    comment : string (default="#")
        Comment character in the begin of a line forces this line to be
        ignored.
    sep : string (default = "\t")
        Specify the separator for the CSV-file.

    Returns
    -------
    l : list
        A list-representation of the CSV file.

    """
    # check for correct fileformat
    if fileformat:
        infile = filename+'.'+fileformat
    else:
        infile = filename
    if not os.path.isfile(infile):
        raise FileNotFoundError(
                "[i] File {0} could not be found.".format(infile)
                )

    l = []
    
    # open the file
    infile = codecs.open(infile,'r','utf-8')

    for line in infile:
        if line.strip() and not line.startswith(comment):
            if strip_lines:
                cells = [c.strip() for c in line.strip().split(sep)]
            else:
                cells = [c for c in line.split(sep)]
            if not dtype:
                l += [cells]
            else:
                l += [[f(c) for f,c in zip(dtype,cells)]]
    infile.close()

    return l

def csv2dict(
        filename,
        fileformat = None,
        dtype = None,
        comment = '#',
        sep = '\t'
        ):
    """
    Very simple function to get quick access to CSV-files.

    Parameters
    ----------
    filename : str
        Name of the input file.
    fileformat : {None str}
        If not specified the file <filename> will be loaded. Otherwise, the
        fileformat is interpreted as the specific extension of the input file.
    dtype : {None list}
        If not specified, all data will be loaded as strings. Otherwise, a
        list specifying the data for each line should be provided.
    comment : string (default="#")
        Comment character in the begin of a line forces this line to be
        ignored.
    sep : string (default = "\t")
        Specify the separator for the CSV-file.

    Returns
    -------
    d : dict
        A dictionary-representation of the CSV file, with the first row being
        used as key and the rest of the rows as values.
    """
    
    l = csv2list(filename,fileformat,dtype,comment,sep)

    d = {}

    for line in l:
        d[line[0]] = line[1:]

    return d               

# define some aliases
def read_csv(
        filename,
        fileformat = '',
        dtype = [],
        comment = '#',
        sep = '\t'
        ):
    """
    Very simple function to get quick access to CSV-files.

    Parameters
    ----------
    filename : str
        Name of the input file.
    fileformat : {None str}
        If not specified the file <filename> will be loaded. Otherwise, the
        fileformat is interpreted as the specific extension of the input file.
    dtype : {None list}
        If not specified, all data will be loaded as strings. Otherwise, a
        list specifying the data for each line should be provided.
    comment : string (default="#")
        Comment character in the begin of a line forces this line to be
        ignored.
    sep : string (default = "\t")
        Specify the separator for the CSV-file.

    Returns
    -------
    l : list
        A list-representation of the CSV file.

    Notes
    -----
    This is but an alias for the csv2dict function.

    """
    return csv2dict(filename,fileformat,dtype,comment,sep)


