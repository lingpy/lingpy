# created: Fr 25 Jan 2013 01:21:36  CET
# modified: Fr 25 Jan 2013 01:21:36  CET

"""
Module provides functions for reading csv-files.
"""

__author__="Johann-Mattis List"
__date__ = "2013-01-25"

def csv2list(
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
    l : list
        A list-representation of the CSV file.

    """

    l = []
    
    if fileformat:
        infile = open(filename+'.'+fileformat)
    else:
        infile = open(filename)
    
    for line in infile:
        if line.strip() and not line.startswith(comment):
            cells = [c.strip() for c in line.strip().split(sep)]
            if not dtype:
                l += [cells]
            else:
                l += [[f(c) for f,c in zip(dtype,cells)]]

    return l

def csv2dict(
        filename,
        fileformat = 'csv',
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
