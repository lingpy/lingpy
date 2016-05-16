# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
import json

from lingpy.algorithm import misc
from lingpy.read.phylip import read_dst, read_scorer
from lingpy.thirdparty import cogent as cg
from lingpy.log import warn, debug
from lingpy.util import read_text_file, setdefaults


def reduce_alignment(alignment):
    """
    Function reduces a given alignment.
    
    Notes
    -----
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
    for i, col in enumerate(cols):
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
        for i, col in enumerate(cols):
            if i not in ignore_indices:
                new_cols += [col]
    else:
        new_cols = cols

    new_alm = misc.transpose(new_cols)

    for i, alm in enumerate(new_alm):
        for j, char in enumerate(alm):
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
        for i, alm in enumerate(alm_clone):
            alm_clone[i] = alm[0].split(' ')
            alm_lens[i] = len(alm_clone[i])

    if len(set(alm_lens)) > 1:
        max_len = max(alm_lens)
        for i, alm in enumerate(alm_clone):
            new_alm = alm + ['-' for x in range(max_len)]
            alm_clone[i] = new_alm[:max_len]

    # then check for alms consisting only of gaps
    cols = misc.transpose(alm_clone)
    idxs = []
    for i, col in enumerate(cols):
        if set(col) == set('-'):
            idxs += [i]
    for idx in idxs[::-1]:
        for i, alm in enumerate(alm_clone):
            del alm_clone[i][idx]
    if alignment != alm_clone:
        lgtxt = 'Modified the alignment:\n'
        for i in range(len(alignment)):
            lgtxt += '[!] ' + ' '.join(alignment[i]) + '->'
            lgtxt += ' '.join(alm_clone[i]) + '\n'
        debug(lgtxt)
        return alm_clone
    else:
        return alignment


def _list2msa(msa_lines, ids=False, header=True, normalize=False, **keywords):
    """
    Function retrieves a dictionary from a list of MSA strings.

    """
    setdefaults(keywords, seq_id='-', dataset='-', input_file='dummy')
    d = dict(ID=[], taxa=[], alignment=[], seqs=[], infile=keywords['input_file'])

    if header:
        start = 2
        d['dataset'] = msa_lines[0]
        d['seq_id'] = msa_lines[1]
    else:
        start = 0
        d['dataset'] = keywords['dataset']
        d['seq_id'] = keywords['seq_id']

    for i, line in enumerate(msa_lines[start:]):
        idx = 1 if ids else 0

        # check for specific id
        if line[0] in ['0', 'LOCAL', 'CROSSED', 'SWAPS', 'MERGE', 'COMPLEX', ]:
            if line[idx] == 'LOCAL':
                d['local'] = []
                for j, x in enumerate(line[idx + 1:]):
                    if x == '*':
                        d['local'] += [j]
            elif line[idx] in ['CROSSED', 'SWAPS']:
                d['swaps'] = []
                swapline = [x for x in line[idx + 1:]]
                j = 0
                while swapline:
                    x = swapline.pop(0)
                    if x == '+':
                        d['swaps'] += [(j, j + 1, j + 2)]
                        swapline.pop(0)
                        swapline.pop(0)
                        j += 2
                    else:
                        pass
                    j += 1
            elif line[idx] in ['COMPLEX', 'MERGE']:
                d['merge'] = {}
                mergeline = [x for x in line[idx + 1:]]
                k = 0
                merge = False
                for j, m in enumerate(mergeline):
                    if m == '<':
                        merge = True
                    if m == '>':
                        merge = False

                    d['merge'][j] = k
                    if not merge:
                        k += 1
            else:
                d[line[idx].lower()] = line[idx + 1:]

        elif line[0] not in ['LOCAL', 'SWAPS', 'MERGE', 'COMPLEX', '0']:
            if ids:
                try:
                    d['ID'] += [int(line[0])]
                except ValueError:
                    d['ID'] += [line[0]]
            else:
                d["ID"] += [i]
            d["taxa"] += [line[idx].rstrip('.')]
            d["seqs"] += [' '.join([l for l in line[idx + 1:] if l != '-'])]
            d["alignment"] += [line[idx + 1:]]

    # normalize the alignment if the option is chosen
    if normalize:
        d['alignment'] = normalize_alignment(d['alignment'])

    return d


def read_msa(infile, comment="#", ids=False, header=True, normalize=True, **keywords):
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

    return _list2msa(msa_lines, header=header, ids=ids, normalize=normalize, **keywords)


def read_qlc(infile, comment='#'):
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
    lines = read_text_file(infile, lines=True, normalize="NFC")
    data, meta, dtype = [], {}, False

    while lines:
        line = lines.pop(0)
        if line.startswith(comment) or not line:
            continue

        if line.startswith('@'):
            key, value = [s.strip() for s in line[1:].split(':', 1)]
            if key == 'tree':
                meta["tree"] = cg.LoadTree(treestring=value)
            elif key == 'json':
                for j1, j2 in json.loads(value).items():
                    meta[j1] = j2
            else:
                if key not in meta:
                    meta[key] = value
                else:
                    if isinstance(meta[key], list):
                        meta[key].append(value)
                    else:
                        warn(
                            "Key '{0}' in input file is not unique! Use JSON-format for "
                            "these datatypes!".format(key))
                        meta[key] = [meta[key]] + [value]
        # line starts with complex stuff
        elif line.startswith('<'):
            tmp = line[1:line.index('>')]
            # check for specific keywords
            if ' ' in tmp:
                dtype = tmp.split(' ')[0]
                keys = {k: v[1:-1]
                        for k, v in [key.split('=') for key in tmp.split(' ')[1:]]}
            else:
                dtype = tmp.strip()
                keys = {}

            tmp = []

            while True:
                line = lines.pop(0)
                if line.startswith('</' + dtype + '>'):
                    break
                tmp += [line]

            tmp = '\n'.join(tmp)

            # check for data stuff
            if dtype == "json":
                tmp = json.loads(tmp)
                if not keys:
                    for key in tmp:
                        meta[key] = tmp[key]
                elif keys:
                    meta[keys["id"]] = {}
                    for k in tmp:
                        meta[keys["id"]][k] = tmp[k]
            elif dtype in ['tre', 'nwk']:
                if "trees" not in meta:
                    meta["trees"] = {}

                if not keys:
                    keys["id"] = "1"

                # XXX consider switching to Tree here XXX
                meta['trees'][keys["id"]] = cg.LoadTree(treestring=tmp)
            elif dtype in ['csv']:
                meta[keys["id"]] = {}
                ncol = int(keys.get('ncol', 2))

                if "dtype" in keys:
                    transf = eval(keys["dtype"])
                else:
                    transf = str

                # split tmp into lines
                tmp = tmp.split('\n')
                for l in tmp:
                    if ncol == 2:
                        a, b = l.split('\t')
                        b = transf(b)
                    else:
                        l = l.split('\t')
                        a = l[0]
                        b = [transf(b) for b in l[1:]]
                    meta[keys["id"]][a] = b
            elif dtype == 'msa':
                tmp = tmp.split('\n')
                if 'msa' not in meta:
                    meta['msa'] = {}

                ref = keys.get('ref', 'cogid')
                if ref not in meta['msa']:
                    meta['msa'][ref] = {}

                tmp_msa = {}
                try:
                    tmp_msa['dataset'] = meta['dataset']
                except:
                    tmp_msa['dataset'] = infile.replace('.csv', '')

                tmp_msa['seq_id'] = keys['id']

                # add consensus string to msa, if it appears in the keys
                if "consensus" in keys:
                    tmp_msa['consensus'] = keys['consensus']

                msad = []
                for l in tmp:
                    if not l.startswith(comment):
                        msad.append([x.strip().rstrip('.') for x in l.split('\t')])
                tmp_msa = _list2msa(msad, header=False, ids=True, **tmp_msa)

                try:
                    meta['msa'][ref][int(keys['id'])] = tmp_msa
                except ValueError:
                    meta['msa'][ref][keys['id']] = tmp_msa

            elif dtype == 'dst':
                taxa, matrix = read_dst(tmp)
                distances = [[0.0 for _ in matrix] for _ in matrix]
                for i, line in enumerate(matrix):
                    for j, cell in enumerate(line):
                        if i < j:
                            distances[i][j] = cell
                            distances[j][i] = cell
                meta['distances'] = distances
            elif dtype == 'scorer':
                scorer = read_scorer(tmp)
                if 'scorer' not in meta:
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
    local_id = data[0][0].lower() in ['id', 'local_id', 'localid']

    # iterate over data and fill the dictionary (a bit inefficient, but enough
    # for the moment)
    try:
        i = 1
        for j, line in enumerate(data[1:]):
            if local_id:
                d[int(line[0])] = line[1:]
            else:
                d[i] = line
                i += 1
    except ValueError as e:
        raise Exception("Error processing line {0}:\n".format(j) +
                        str(data[1:][j]) + '\nOriginal error message: ' + str(e))

    # assign the header to d[0]
    if local_id:
        d[0] = [x.lower() for x in data[0][1:]]
    else:
        d[0] = [x.lower() for x in data[0]]

    for m in meta:
        d[m] = meta[m]

    if 'trees' in d and 'tree' not in d:
        d['tree'] = sorted(d['trees'].items(), key=lambda x: x[0])[0][1]

    return d
