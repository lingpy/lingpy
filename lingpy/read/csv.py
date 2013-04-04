# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 17:02
# modified : 2013-03-05 08:37
"""
Module provides functions for reading csv-files.
"""

__author__="Johann-Mattis List"
__date__="2013-03-05"

import json

# lingpy-internal imports
from ..thirdparty import cogent as cg
from .phylip import read_dst

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

def qlc2dict_deprecated(infile):
    """
    Simple function that loads qlc-format into a dictionary.

    Parameters
    ----------
    infile : str
        Name of the input file.

    Returns
    -------
    d : dict
        A dictionary with integer keys corresponding to the order of the lines
        of the input file. The header is given 0 as a specific key.
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
    try:
        i = 1
        for line in data[1:]:
            if local_id:
                d[int(line[0])] = line[1:]
            else:
                d[i] = line
                i += 1
    except:
        print("[!] Something is wrong with your input file. If it contains an ID column, make sure it consists only of integers.")

    # assign the header to d[0]
    if local_id:
        d[0] = [x.lower() for x in data[0][1:]]
    else:
        d[0] = [x.lower() for x in data[0]]

    # return the stuff
    return d

# define some aliases
def read_csv(
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

    Notes
    -----
    This is but an alias for the csv2dict function.

    """
    return csv2dict(filename,fileformat,dtype,comment,sep)

def read_qlcOld(
        infile
        ):
    """
    Simple function that loads qlc-format into a dictionary.

    Parameters
    ----------
    infile : str
        Name of the input file.

    Returns
    -------
    d : dict
        A dictionary with integer keys corresponding to the order of the lines
        of the input file. The header is given 0 as a specific key.
    
    Notes
    -----
    This is but an alias for the qlc2dict function.
    """
    return qlc2dict(infile)

def read_qlc(infile):
    """
    Simple function that loads qlc-format into a dictionary.

    Parameters
    ----------
    infile : str
        Name of the input file.

    Returns
    -------
    d : dict
        A dictionary with integer keys corresponding to the order of the lines
        of the input file. The header is given 0 as a specific key.
    """
    
    # open the file
    try:
        inf = open(infile)
    except:
        raise ValueError("[i] Infile could not be opened.")
    
    # create data array
    data = []
    meta = {}

    # read lines from infile
    lines = inf.readlines()

    # set dtype to false
    dtype = False

    # start the loop
    while lines:
        line = lines.pop(0)

        # line is empty or starts with hash
        if line.startswith('#') or not line.strip():
            pass
        
        # line starts with key-value @
        elif line.startswith('@'):
            key = line[1:line.index(':')]
            value = line[line.index(':')+1:].strip()
            if key not in ['tree']:
                meta[key] = value
            else:
                if key == 'tree':
                    meta["tree"] = cg.LoadTree(treestring=value)
        
        # line starts with complex stuff
        elif line.startswith('<'):
            tmp = line[1:line.index('>')]
            # check for specific keywords
            if ' ' in tmp:
                dtype = tmp.split(' ')[0]
                keys = dict([(k,v[1:-1]) for k,v in [key.split('=') for key in
                    tmp.split(' ')[1:]]])
            else:
                dtype = tmp.strip()
                keys = {}

            tmp = []

            # iterate
            while True:
                line = lines.pop(0)

                if line.startswith('</'+dtype+'>'):
                    break
                else:
                    tmp += [line.strip()]

            tmp = '\n'.join(tmp)

            # check for data stuff
            # json
            if dtype == "json":
                tmp = json.loads(tmp)
                if not keys:
                    for key in tmp:
                        meta[key] = tmp[key]
                elif keys:
                    meta[keys["id"]] = {}
                    for k in tmp:
                        meta[keys["id"]][k] = tmp[k]
            
            # tree
            elif dtype in ['tre','nwk']:
                # check for "tree" in meta
                if "tree" not in meta:
                    meta["tree"] = {}
                
                # add the data
                if not keys:
                    keys["id"] = "t1"

                meta['tree'][keys["id"]] = cg.LoadTree(treestring=tmp)

            # csv
            elif dtype in ['csv']:

                meta[keys["id"]] = {}

                # check for columns
                if "ncol" in keys:
                    ncol = int(keys["ncol"])
                else:
                    ncol = 2

                # check for dtype
                if "dtype" in keys:
                    transf = eval(keys["dtype"])
                else:
                    transf = str

                # split tmp into lines
                tmp = tmp.split('\n')
                for l in tmp:
                    if ncol == 2:
                        a,b = l.split('\t')
                        b = transf(b)
                    else:
                        l = l.split('\t')
                        a = l[0]
                        b = [transf(b) for b in l[1:]]
                    meta[keys["id"]][a] = b
            elif dtype == 'msa':
                
                # check for dtype
                tmp = tmp.split('\n')
                if 'msa' not in meta:
                    meta['msa'] = {}
                tmp_msa = {}
                try:
                    tmp_msa['dataset'] =  meta['dataset']
                except:
                    tmp_msa['dataset'] = infile.replace('.csv','')
                tmp_msa['seq_id'] = keys['id']
                tmp_msa['seqs'] = []
                tmp_msa['alignment'] = []
                tmp_msa['taxa'] = []
                for l in tmp[2:]:
                    this_line = l.split('\t')
                    tmp_msa['taxa'] += [this_line[0]]
                    tmp_msa['seqs'] += [' '.join(this_line[1:])]
                    tmp_msa['alignment'] += [this_line[1:]]
                try:
                    meta['msa'][int(keys['id'])] = tmp_msa   
                except:
                    meta['msa'][keys['id']] = tmp_msa
            
            elif dtype == 'dst':
                taxa,matrix = read_dst(tmp)
                new_taxa = sorted(taxa)
                distances = [[0.0 for line in matrix] for line in matrix]
                for i,line in enumerate(matrix):
                    for j,cell in enumerate(line):
                        if i < j:
                            idxA = new_taxa.index(taxa[i])
                            idxB = new_taxa.index(taxa[j])
                            distances[i][j] = cell
                            distances[j][i] = cell
                meta['distances'] = distances
            elif dtype == 'taxa':
                meta['taxa'] = [t.strip() for t in tmp.split('\n')]          
        else:
            data += [[l.strip() for l in line.split('\t')]]
    
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
    try:
        i = 1
        for line in data[1:]:
            if local_id:
                d[int(line[0])] = line[1:]
            else:
                d[i] = line
                i += 1
    except:
        print("[!] Something is wrong with your input file. If it contains an ID column, make sure it consists only of integers.")

    # assign the header to d[0]
    if local_id:
        d[0] = [x.lower() for x in data[0][1:]]
    else:
        d[0] = [x.lower() for x in data[0]]
    
    for m in meta:
        d[m] = meta[m]

    return d

def qlc2dict(infile):

    return read_qlc(infile)



