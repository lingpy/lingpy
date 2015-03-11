# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

from ..algorithm import misc
from .csv import csv2list
from .phylip import read_dst,read_scorer
from ..thirdparty import cogent as cg
from ..log import warn
import json
import os
import unicodedata
from ..util import read_text_file

def reduce_alignment(alignment):
    """
    Function reduces a given alignment.
    
    Note
    ----
    Reduction here means that the output alignment consists only of those parts
    which have not been marked to be ignored by the user (parts in brackets).
    It requires that all data is properly coded. If reduction fails, this will
    throw a warning, and all brackets are simply removed in the output
    alignment.
    """
    
    # check for bracket indices in all columns
    cols = misc.transpose(alignment)
    
    ignore_indices = []
    ignore = False
    for i,col in enumerate(cols):
        reduced_col = sorted(set(col))
        
        if '(' in reduced_col:
            if len(reduced_col) == 1:
                ignore_indices += [i] 
                ignore = True
            else:
                ignore = False
        elif ')' in reduced_col:
            if len(reduced_col) == 1:
                ignore_indices += [i]
                ignore = False
            else:
                ignore_indices = []
        elif ignore:
            ignore_indices += [i]

    if ignore_indices:
        new_cols = []
        for i,col in enumerate(cols):
            if i not in ignore_indices:
                new_cols += [col]
    else:
        new_cols = cols

    new_alm = misc.transpose(new_cols)

    for i,alm in enumerate(new_alm):
        for j,char in enumerate(alm):
            if char in '()':
                new_alm[i][j] = '-'

    return new_alm

def normalize_alignment(alignment):
    """
    Function normalizes an alignment.

    Normalization here means that columns consisting only of gaps will be
    deleted, and all sequences will be stretched to equal length by adding
    additional gap characters in the end of smaller sequences.
    """
    # clone the alignment
    alm_clone = [[x for x in y] for y in alignment]

    # first check for alms of different length
    alm_lens = [len(alm) for alm in alm_clone]
    if alm_lens.count(1) == len(alm_lens):
        for i,alm in enumerate(alm_clone):
            alm_clone[i] = alm[0].split(' ')
            alm_lens[i] = len(alm_clone[i])

    if len(set(alm_lens)) > 1:
        max_len = max(alm_lens)
        for i,alm in enumerate(alm_clone):
            new_alm = alm + ['-' for x in range(max_len)]
            alm_clone[i] = new_alm[:max_len]

    # then check for alms consisting only of gaps
    cols = misc.transpose(alm_clone)
    idxs = []
    for i,col in enumerate(cols):
        if set(col) == set('-'):
            idxs += [i]
    for idx in idxs[::-1]:
        for i,alm in enumerate(alm_clone):
            del alm_clone[i][idx]
    if alignment != alm_clone:
        lgtxt = '[!] Modified the alignment:\n'
        for i in range(len(alignment)):
            lgtxt += '[!] '+ ' '.join(alignment[i]) + '->'
            lgtxt += ' '.join(alm_clone[i]) + '\n'
        warn(lgtxt)
        return alm_clone
    else:
        return alignment

def _list2msa(
        msa_lines,
        ids=False,
        header=True,
        normalize = False,
        **keywords
        ):
    """
    Function retrieves a dictionary from a list of MSA strings.

    """
    defaults = dict(
            seq_id = '-',
            dataset = '-',
            input_file = 'dummy'
            )
    for key in defaults:
        if key not in keywords:
            keywords[key] = defaults[key]

    # create the dictionary
    d = dict(
            ID = [],
            taxa = [],
            alignment = [],
            seqs = [],
            infile = keywords['input_file']
            )
    
    if header:
        start = 2
        d['dataset'] = msa_lines[0]
        d['seq_id'] = msa_lines[1]
    else:
        start = 0
        d['dataset'] = keywords['dataset']
        d['seq_id'] = keywords['seq_id']

    for i,line in enumerate(msa_lines[start:]):
        if ids:
            idx = 1
        else:
            idx = 0
        
        # check for specific id
        if line[0] in  ['0', 'LOCAL', 'CROSSED', 'SWAPS', 'MERGE', 'COMPLEX',]:
            if line[idx] == 'LOCAL':
                d['local'] = []
                for j,x in enumerate(line[idx+1:]):
                    if x == '*':
                        d['local'] += [j]
            elif line[idx] in ['CROSSED', 'SWAPS']:
                d['swaps'] = []
                swapline = [x for x in line[idx+1:]]
                j=0
                while swapline:
                    x = swapline.pop(0)
                    if x == '+':
                        d['swaps'] += [(j,j+1,j+2)]
                        swapline.pop(0)
                        swapline.pop(0)
                        j += 2
                    else:
                        pass
                    j += 1
            elif line[idx] in ['COMPLEX', 'MERGE']:
                d['merge'] = {}
                mergeline = [x for x in line[idx+1:]]
                k = 0
                merge = False
                for j,m in enumerate(mergeline):

                    if m == '<':
                        merge = True
                    if m == '>':
                        merge = False

                    d['merge'][j] = k
                    if not merge:
                        k += 1

            else:
                d[line[idx].lower()] = line[idx+1:]

        elif line[0] not in ['LOCAL','SWAPS','MERGE','COMPLEX','0']:
            if ids:
                try:
                    d['ID'] += [int(line[0])]
                except:
                    d['ID'] += [line[0]]
            else:
                d["ID"] += [i]
            d["taxa"] += [line[idx].rstrip('.')]
            d["seqs"] += [' '.join([l for l in line[idx+1:] if l != '-'])]
            d["alignment"] += [line[idx+1:]]
    
    # normalize the alignment if the option is chosen
    if normalize:
        d['alignment'] = normalize_alignment(d['alignment'])

    return d

def read_msa(
        infile,
        comment = "#",
        ids = False,
        header = True,
        normalize = True,
        **keywords
        ):
    """
    Simple function to load an MSA object.

    Parameters
    ----------
    infile : str
        The name of the input file.
    comment : str (default="#")
        The comment character. If a line starts with this character, it will be
        ignored.
    ids : bool (default=False)
        Indicate whether the MSA file contains unique IDs for all sequences or
        not.

    Returns
    -------
    d : dict
        A dictionary in which keys correspond to specific parts of a multiple
        alignment. This dictionary can be directly passed to alignment
        functions, such as :py:class:`lingpy.sca.MSA`.
    """
    if 'input_file' not in keywords:
        keywords['input_file'] = infile

    f = read_text_file(infile, normalize='NFC', lines=True)
    msa_lines = []
    for line in f:
        if line.strip() and not line.startswith(comment):
            newlines = [t.strip().rstrip('.') for t in line.split('\t')]
            if len(newlines) == 1:
                msa_lines += newlines
            else:
                msa_lines += [newlines]
    
    d = _list2msa(msa_lines, header=header, ids=ids, normalize=normalize,
            **keywords)

    return d
    
def read_qlc(
        infile,
        comment = '#'
        ):
    """
    Simple function that loads qlc-format into a dictionary.

    Parameters
    ----------
    infile : str
        The name of the input file.
    comment : str (default="#")
        The comment character. If a line starts with this character, it will be
        ignored.

    Returns
    -------
    d : dict
        A dictionary with integer keys corresponding to the order of the lines
        of the input file. The header is given 0 as a specific key.
    """
    # read lines from text file
    lines = read_text_file(infile, lines=True, normalize="NFC")
    
    # create data array
    data = []
    meta = {}

    # set dtype to false
    dtype = False

    # start the loop
    while lines:
        line = lines.pop(0)

        # line is empty or starts with hash
        if line.startswith(comment) or not line:
            pass
        
        # line starts with key-value @
        elif line.startswith('@'):
            key = line[1:line.index(':')]
            value = line[line.index(':')+1:].strip()
            if key not in ['tree','json']:
                if key not in meta:
                    meta[key] = value
                else:
                    if type(meta[key]) == list:
                        meta[key] += [value]
                    else:
                        print("[WARNING] Key '{0}' in input file is not unique! Use JSON-format for these datatypes!".format(key))
                        meta[key] = [meta[key]]
                        meta[key] += [value]
                        
            else:
                if key == 'tree':
                    meta["tree"] = cg.LoadTree(treestring=value)
                if key == 'json':
                    jsonpairs = json.loads(value)
                    for j1,j2 in jsonpairs.items():
                        meta[j1] = j2
        
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
                    tmp += [line]

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
                if "trees" not in meta:
                    meta["trees"] = {}
                
                # add the data
                if not keys:
                    keys["id"] = "1"
                
                # XXX consider switching to Tree here XXX
                meta['trees'][keys["id"]] = cg.LoadTree(treestring=tmp)

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

                # check for additional reference attribute
                if 'ref' in keys:
                    ref = keys['ref']
                else:
                    ref = 'cogid' # XXX check this later for flexibility

                # check for msa.ref:
                if ref not in meta['msa']:
                    meta['msa'][ref] = {}

                tmp_msa = {}
                try:
                    tmp_msa['dataset'] =  meta['dataset']
                except:
                    tmp_msa['dataset'] = infile.replace('.csv','')

                tmp_msa['seq_id'] = keys['id']

                # add consensus string to msa, if it appears in the keys
                if "consensus" in keys:
                    tmp_msa['consensus'] = keys['consensus']
                
                msad = []
                for l in tmp:
                    if not l.startswith(comment):
                        
                        line = [x.strip().rstrip('.') for x in
                                l.split('\t')]
                        msad += [line]
                tmp_msa = _list2msa(msad, header=False, ids=True, **tmp_msa)

                try:
                    meta['msa'][ref][int(keys['id'])] = tmp_msa   
                except:
                    meta['msa'][ref][keys['id']] = tmp_msa
            
            elif dtype == 'dst':
                taxa, matrix = read_dst(tmp)
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
            elif dtype == 'scorer':
                scorer = read_scorer(tmp)
                if not 'scorer' in meta:
                    meta['scorer'] = {}
                if 'id' not in keys:
                    keys['id'] = 'basic'
                meta['scorer'][keys['id']] = scorer

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
        for j,line in enumerate(data[1:]):
            if local_id:
                d[int(line[0])] = line[1:]
            else:
                d[i] = line
                i += 1
    except ValueError as e:
        raise Exception("Error processing line {0}:\n".format(j) +
                str(data[1:][j])+'\nOriginal error message: '+str(e))

    # assign the header to d[0]
    if local_id:
        d[0] = [x.lower() for x in data[0][1:]]
    else:
        d[0] = [x.lower() for x in data[0]]
    
    for m in meta:
        d[m] = meta[m]
    
    # check for tree-attributes
    if 'trees' in d and not 'tree' in d:
        d['tree'] = sorted(
                d['trees'].items(),
                key = lambda x: x[0],
                )[0][1]
                

    return d
