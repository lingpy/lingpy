"""
Module provides functions for reading csv-files.
"""
from __future__ import unicode_literals, division, print_function

from lingpy.util import read_text_file
from lingpy.sequence.sound_classes import asjp2tokens


def csv2list(
    filename,
    fileformat='',
    dtype=None,
    comment='#',
    sep='\t',
    strip_lines=True,
    header=False
):
    """
    Very simple function to get quick (and somewhat naive) access to CSV-files.

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
        infile = filename + '.' + fileformat
    else:
        infile = filename

    if dtype is None:
        dtype = []

    l = []

    # open the file
    infile = read_text_file(infile, lines=True, normalize="NFC")

    # check for header
    idx = 0 if header else -1

    for i, line in enumerate(infile):
        if line and not line.startswith(comment) and idx != i:
            if strip_lines:
                cells = [c.strip() for c in line.strip().split(sep)]
            else:
                cells = [c.strip() for c in line.split(sep)]
            if not dtype:
                l += [cells]
            else:
                l += [[f(c) for f, c in zip(dtype, cells)]]

    return l


def csv2dict(
    filename,
    fileformat=None,
    dtype=None,
    comment='#',
    sep='\t',
    strip_lines=True,
    header=False
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
    d : dict
        A dictionary-representation of the CSV file, with the first row being
        used as key and the rest of the rows as values.
    """

    l = csv2list(filename, fileformat, dtype, comment, sep, strip_lines, header)
    return {line[0]: line[1:] for line in l}


def csv2multidict(filename, comment='#', sep='\t'):
    """
    Function reads a csv-file into a multi-dimensional dictionary structure.
    """
    tsv = csv2list(filename, comment=comment, sep=sep)
    header = tsv[0]
    return {line[0]: dict(zip(header[1:], line[1:])) for line in tsv[1:]}


def read_asjp(
    infile,
    family='Indo-European',
    classification='hh',
    max_synonyms=2,
    min_population=lambda x: x > 0 or abs(x) > 1900,
    merge_vowels=True,
    evaluate=False
):
    # read in all data
    data = csv2list(infile, strip_lines=False)

    # find the data confirming to selection, get the header for the
    # classification first
    header = [h.lower() for h in data[0]]

    # check for family type
    if not evaluate:
        evaluate = lambda x, y, z: x[y].startswith(z)

    # index the classification index
    if ',' in classification:
        a, b = classification.split(',')
        clsA_idx = header.index(a)
        clsB_idx = header.index(b)
        cls_idx = (clsA_idx, clsB_idx)
    else:
        cls_idx = header.index(classification)

    # create dictionary to store the data
    D, idx = {}, 1

    # lats / longs
    meta = dict(
        coords={},
        iso={},
        population={},
        classification=dict(
            hammarstroem={},
            ethnologue={},
            wals={},
            wals_genus={}
        )
    )

    # iterate over data and extract the lines
    for line in data[1:]:
        if evaluate(line, cls_idx, family):
            lang = line[0].strip()
            wls = line[1].strip()
            wls_gen = line[2].strip()
            eth = line[3].strip()
            hh = line[4].strip()
            try:
                lat = float(line[5].replace(',', '.'))
                lng = float(line[6].replace(',', '.'))
            except:
                lat = ''
                lng = ''
            pop = int(line[7]) if line[7] else -10
            iso = line[9].strip()

            # check for population
            if min_population(pop):  # >= min_population:
                # append data to meta
                if lang not in meta['coords']:
                    meta["coords"][lang] = (lat, lng)
                    meta["classification"]["hammarstroem"][lang] = hh
                    meta["classification"]["wals"][lang] = wls
                    meta["classification"]["wals_genus"][lang] = wls_gen
                    meta["classification"]["ethnologue"][lang] = eth
                    meta["population"][lang] = pop
                    meta["iso"][lang] = iso

                    for i, items in enumerate(line[10:], 10):
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
                                entry = entry.replace(' ', '_')
                            tokens = ' '.join(asjp2tokens(entry))
                            D[idx] = [lang, item, entry, tokens.split(' '), loan]
                            idx += 1
    D[0] = ['doculect', 'concept', 'counterpart', 'tokens', 'known_borrowings']
    D['meta'] = meta
    return D
