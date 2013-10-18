# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 17:02
# modified : 2013-10-18 23:31
"""
Module provides functions for reading csv-files.
"""

__author__="Johann-Mattis List"
__date__="2013-10-18"

import codecs
import os
from ..settings import rcParams
from ..sequence.sound_classes import ipa2tokens
import re

def csv2list(
        filename,
        fileformat = '',
        dtype = [],
        comment = '#',
        sep = '\t',
        strip_lines = True,
        header = False
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
    strip_lines : bool (default=True)
        Specify whether empty "cells" in the input file should be preserved. If
        set to c{False}, each line will be stripped first, and all whitespace
        will be cleaned. Otherwise, each line will be separated using the
        specified separator, and no stripping of whitespace will be carried
        out.
    header : bool (default=False)
        Indicate, whether the data comes along with a header.

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
        raise NameError(
                "[ERROR] File {0} could not be found.".format(infile)
                )

    l = []
    
    # open the file
    infile = codecs.open(infile,'r','utf-8')
    
    # check for header
    if header:
        idx = 0
    else:
        idx = -1

    for i,line in enumerate(infile):
        if line.strip() and not line.startswith(comment) and idx != i:
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
        sep = '\t',
        strip_lines = True,
        header = False
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
    header : bool (default=False)
        Indicate, whether the data comes along with a header.

    Returns
    -------
    d : dict
        A dictionary-representation of the CSV file, with the first row being
        used as key and the rest of the rows as values.
    """
    
    l = csv2list(filename,fileformat,dtype,comment,sep,strip_lines,header)

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

def read_asjp(
        infile,
        family = 'Indo-European',
        classification = 'hh',
        max_synonyms = 2,
        min_population = 0,
        merge_vowels=True,
        evaluate = False
        ):
    
    # read in all data
    data = csv2list(infile,strip_lines=False)
    
    # find the data confirming to selection, get the header for the
    # classification first
    header = [h.lower() for h in data[0]]
    
    # check for family type
    if not evaluate:
        evaluate = lambda x:x.startswith(family)

    # index the classification index
    cls_idx = header.index(classification)
    
    # create dictionary to store the data
    D,idx = {},1

    # lats / longs
    meta = dict(
            coords = {},
            iso = {},
            population = {},
            classification = dict(
                hammarstroem = {},
                ethnologue = {},
                wals = {},
                wals_genus = {}
                )
            )

    # iterate over data and extract the lines
    for line in data[1:]:
        if evaluate(line[cls_idx]):
            lang = line[0].strip()
            wls = line[1].strip()
            wls_gen = line[2].strip()
            eth = line[3].strip()
            hh = line[4].strip()
            try:
                lat = float(line[5].replace(',','.'))
                lng = float(line[6].replace(',','.'))
            except:
                lat = ''
                lng = ''
            pop = int(line[7]) if line[7] else ''
            iso = line[9].strip()
            
            # check for population
            if pop >= min_population:
                # append data to meta
                meta["coords"][lang] = (lat,lng)
                meta["classification"]["hammarstroem"][lang] = hh
                meta["classification"]["wals"][lang] =  wls
                meta["classification"]["wals_genus"][lang] = wls_gen
                meta["classification"]["ethnologue"][lang] = eth
                meta["population"][lang] = pop
                meta["iso"][lang] = iso

                for i,items in enumerate(line[10:],10):

                    item = header[i].strip()
                    entries = [e.strip() for e in line[i].split(',') if
                            e.strip() and 'xxx' not in e.lower()][:max_synonyms]
                
                    for entry in entries:
                        if entry.startswith('%'):
                            entry = entry[1:]
                            loan = 1
                        else:
                            loan = 0
                        if ' ' in entry:
                            entry = entry.replace(' ','_')
                        tokens = ' '.join(
                                    ipa2tokens(
                                        entry,
                                        diacritics = '*$~',
                                        vowels = 'aeiouE3',
                                        tones = '',
                                        combiners = '',
                                        merge_vowels = merge_vowels#rcParams['merge_vowels']
                                        )
                                    )
                        tokens = re.sub(r'([^ ]) ([^ ])~',r'\1\2~',tokens)

                        D[idx] = [
                                lang,
                                item,
                                entry,
                                tokens.split(' '),
                                loan
                                ]
                        idx += 1
    D[0] = ['doculect','concept','counterpart','tokens','known_borrowings']
    D['meta'] = meta

    return D
                                

                    



