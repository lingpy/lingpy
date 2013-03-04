# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 17:02
# modified : 2013-03-04 17:02
"""
Module provides functions for reading csv-files.
"""

__author__="Johann-Mattis List"
__date__="2013-03-04"

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

def qlc2dict(infile):
    """
    Simple function that loads qlc-format into a dictionary.

    """

    # read the data into a list
    data = []

    # open the file
    f = open(infile)

    for line in f:
        # ignore hashed lines
        if not line.startswith('#') and not line.startswith('@'):

            # mind to strip newlines
            data.append(line.strip('\n\r').split('\t'))
    
    # create the dictionary in which the data will be stored
    d = {}

    # check for first line, if a local ID is given in the header (or simply
    # "ID"), take this line as the ID, otherwise create it
    if data[0][0].lower() in ['id','local_id','localid']:
        local_id = True
    else:
        local_id = False

    # iterate over data and fill the dictionary (a bit inefficient, but enough
    # for the moment)
    i = 1
    for line in data[1:]:
        if local_id:
            d[int(line[0])] = line[1:]
        else:
            d[i] = line
            i += 1

    # assign the header to d[0]
    if local_id:
        d[0] = [x.lower() for x in data[0][1:]]
    else:
        d[0] = [x.lower() for x in data[0]]

    # return the stuff
    return d

