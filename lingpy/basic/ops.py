"""
Module provides basic operations on Wordlist-Objects.
"""
from __future__ import unicode_literals, print_function, absolute_import, division

import json
from string import ascii_letters, digits
from collections import defaultdict
from itertools import product

from six import text_type

from lingpy.settings import rcParams
from lingpy.convert.strings import matrix2dst, scorer2str, msa2str
from lingpy.algorithm import clustering, misc
from lingpy import util
from lingpy import log


def get_score(wl, ref, mode, taxA, taxB, concepts_attr='concepts',
        ignore_missing=False):
    if mode in ['shared', 'jaccard']:
        listA, listB = [wl.get_list(col=tax, entry=ref) for tax in [taxA, taxB]]
        shared = [x for x in listA if x in listB]

        if mode == 'jaccard':
            return 1 - len(set(shared)) / len(set(listA + listB))
        return len(shared)

    assert mode == 'swadesh'
    # get the two dictionaries
    dictA, dictB = [wl.get_dict(col=tax, entry=ref) for tax in [taxA, taxB]]

    # count amount of shared concepts
    shared, missing = 0, 0

    for concept in getattr(wl, concepts_attr):
        if concept not in dictA or concept not in dictB:
            missing += 1 if not ignore_missing else 0
        elif [k for k in dictA[concept] if k in dictB[concept]]:
            shared += 1

    try:
        return 1 - shared / (wl.height - missing)
    except ZeroDivisionError:
        log.get_logger().exception(
            "Zero-division error encountered in '{0}' and '{1}'.".format(taxA, taxB))
        return 1.0


def wl2dst(
        wl,  # wordlist object
        taxa="taxa",
        concepts="concepts",
        ref='cogid',
        refB='',
        mode='swadesh',
        ignore_missing=False,
        **keywords):
    """
    Function converts wordlist to distance matrix.
    """
    # check for attributes
    assert hasattr(wl, taxa) and hasattr(wl, concepts)

    distances = [[0 for i in range(wl.width)] for j in range(wl.width)]

    for (i, taxA), (j, taxB) in product(enumerate(getattr(wl, taxa)), repeat=2):
        if i < j:
            score = get_score(wl, ref, mode, taxA, taxB,
                    concepts_attr=concepts, ignore_missing=ignore_missing)
            distances[i][j] = score
            if not refB:
                distances[j][i] = score
        elif i == j:
            if mode == 'shared':
                distances[i][j] = len(wl.get_list(col=taxA, flat=True))
        elif i > j and refB:
            distances[i][j] = get_score(
                wl, refB, mode, taxA, taxB, concepts_attr=concepts,
                ignore_missing=ignore_missing)

    return distances


def wl2dict(
        wordlist,
        sections,
        entries,
        exclude=None):
    """
    Convert a wordlist to a complex dictionary with headings as keys.
    """
    assert sections

    # define output dictionary
    out = {}
    exclude = exclude or []

    # determine the last section
    sorted_sections = sorted(sections)
    last_section = sorted_sections[-1]

    for key in wordlist:
        if key not in exclude:
            # pass temporary pointer
            tmp = out

            for s in sorted_sections:
                # get datapoint and text
                data_point = wordlist[key, sections[s][0]]
                outkey = (data_point, sections[s][1].format(data_point))

                # dive deeper if this is not the last section
                if s != last_section:
                    # dive deeper
                    if outkey not in tmp:
                        tmp[outkey] = {}
                    tmp = tmp[outkey]
                else:
                    # dive to last level
                    if outkey not in tmp:
                        tmp[outkey] = []
                    tmp = tmp[outkey]

                    # get the final list of entries
                    tmp_list = []
                    for entry, format_string in entries:
                        if type(entry) in (list, tuple):
                            entry = ' '.join(entry)
                        tmp_list.append(format_string.format(wordlist[key, entry]))
                    tmp += [tmp_list]

    return out


def renumber(wordlist, source, target='', override=False):
    """
    Create numerical identifiers from string identifiers.
    """
    # iterate over wordlist and get all source ids
    sources = sorted(set([text_type(wordlist[k, source]) for k in wordlist]))

    # convert to numbers
    targets = list(range(1, len(sources) + 1))

    # add to wordlist
    target = target or (source + 'id')

    # make converter
    converter = dict(zip(sources, targets))

    # check for zero ids
    if 0 in converter:
        converter[0] = 0
    if '' in converter:
        converter[''] = 0

    wordlist.add_entries(
        target, source, lambda x: converter[text_type(x)], override=override)

    # add stuff to meta
    wordlist._meta[source + '2' + target] = converter
    log.info("Successfully renumbered {0}.".format(source))


def clean_taxnames(
        wordlist,
        column='doculect',
        f=lambda x: ''.join([t for t in x if t not in '()[]{},;:'])
        .replace('-', '_').replace(' ', '_')):
    """
    Function cleans taxon names in order to make sure they can be used in Newick files.

    """
    # clean the names for all taxa in a wordlist
    current_taxa = eval('wordlist.' + column)
    new_taxa = [f(taxon) for taxon in current_taxa]

    old2new = dict(zip(current_taxa, new_taxa))
    new2old = dict(zip(new_taxa, current_taxa))

    if column == wordlist._col_name:
        wordlist.cols = [old2new[t] for t in current_taxa]

    wordlist.add_entries('_doculect', 'doculect', lambda x: old2new[x], override=True)
    wordlist.add_entries('doculect', '_doculect', lambda x: new2old[x], override=True)

def calculate_data(
        wordlist,
        data,
        taxa='taxa',
        concepts='concepts',
        ref='cogid',
        **keywords):
    """
    Manipulate a wordlist object by adding different kinds of data.

    Parameters
    ----------
    data : str
        The type of data that shall be calculated. Currently supports

        * "tree": calculate a reference tree based on shared cognates
        * "dst": get distances between taxa based on shared cognates
        * "cluster": cluster the taxa into groups using different methods


    """
    logger = log.get_logger()
    util.setdefaults(
        keywords,
        distances=False,
        tree_calc="upgma",
        cluster="upgma",
        force=False,
        threshold=0.5,
        cluster_method='upgma')

    # get taxa for current calculation
    these_taxa = eval('wordlist.' + taxa)

    # calculate distances
    if data in ['distances', 'dst']:
        wordlist._meta['distances'] = wl2dst(wordlist, taxa, concepts, ref, **keywords)
    elif data in ['diversity', 'div']:
        etd = wordlist.get_etymdict(ref=ref)
        wordlist._meta['diversity'] = \
            (len(etd) - wordlist.height) / (len(wordlist) - wordlist.height)
    elif data in ['tre', 'tree', 'nwk']:
        if 'distances' not in wordlist._meta:
            wordlist._meta['distances'] = \
                wl2dst(wordlist, taxa, concepts, ref, **keywords)
        distances = wordlist._meta['distances']
        if 'tree' in wordlist._meta and not keywords['force']:
            logger.warn("Reference tree has already been calculated, force overwrite by "
                        "setting 'force' to 'True'.")
            return
        wordlist._meta['tree'] = clustering.matrix2tree(
            distances, these_taxa, keywords['tree_calc'], keywords['distances'])

    elif data in ['groups', 'cluster']:
        if 'distances' not in wordlist._meta:
            distances = wl2dst(wordlist, taxa, concepts, ref, **keywords)
        else:
            distances = wordlist._meta['distances']
        if 'groups' in wordlist._meta and not keywords['force']:
            logger.warn("Distance matrix has already been calculated, force overwrite by "
                        "setting 'force' to 'True'.")
            return
        wordlist._meta['groups'] = clustering.matrix2groups(
            keywords['threshold'], distances, these_taxa, keywords['cluster_method'])
    log.info("Successfully calculated {0}.".format(data))


def wl2qlc(
        header,
        data,
        filename='',
        formatter='concept',
        **keywords):
    """
    Write the basic data of a wordlist to file.
    """
    util.setdefaults(
        keywords,
        ignore=['taxa', 'doculects', 'msa'],
        fileformat='qlc',
        prettify=True)
    if keywords['ignore'] == 'all':
        keywords['ignore'] = [
            'taxa', 'scorer', 'meta', 'distances', 'doculects', 'msa', 'json']

    formatter = formatter.upper()
    if not filename:
        filename = rcParams['filename']

    # create output string
    out = '# Wordlist\n' if keywords['prettify'] else ''

    # write meta to file
    meta = keywords.get("meta", {})
    kvpairs = {}
    jsonpairs = {}
    msapairs = {}
    trees = {}
    distances = ''
    taxa = ''
    scorer = ''

    for k, v in meta.items():
        # simple key-value-pairs
        if isinstance(v, (text_type, int)) or k == "tree":
            kvpairs[k] = v
        elif k == 'msa' and k not in keywords['ignore']:
            # go a level deeper, checking for keys
            for ref in v:
                if ref not in msapairs:
                    msapairs[ref] = {}
                for a, b in v[ref].items():
                    msapairs[ref][a] = b
        elif k == 'distances':
            distances = matrix2dst(v, meta['taxa'])
        elif k in ['taxa', 'doculect', 'taxon', 'doculects']:
            # we need to find a better solution here, since it is not nice to
            # have taxa written to json again and again
            pass
        elif k == 'trees' and k not in keywords['ignore']:
            trees = ''
            for key, value in v.items():
                trees += '<tre id="{0}">\n{1}\n</tre>\n'.format(key, value)
        elif k == 'scorer' and k not in keywords['ignore']:
            for key, value in v.items():
                scorer += '<{2} id="{0}">\n{1}</{2}>\n\n'.format(
                    key, scorer2str(value), k)
        else:
            # check whether serialization works
            try:
                json.dumps(v)
                jsonpairs[k] = v
            except TypeError:
                pass

    if kvpairs and 'meta' not in keywords['ignore']:
        out += '\n# META\n' if keywords['prettify'] else ''
        for k, v in sorted(kvpairs.items(), key=lambda x: x[0]):
            out += '@{0}:{1}\n'.format(k, v)
    if taxa and keywords['taxa']:
        out += '\n# TAXA\n<taxa>\n' + taxa + '\n</taxa>\n'
    if jsonpairs and 'json' not in keywords['ignore']:
        out += "@json: " + json.dumps(jsonpairs) + '\n'
    if msapairs and 'msa' not in keywords['ignore']:
        for ref in msapairs:
            out += "\n# MSA reference: {0}\n".format(ref)
            for k, v in msapairs[ref].items():
                if 'consensus' in v:
                    out += '#\n<msa id="{0}" ref="{1}" consensus="{2}">\n'.format(
                        k, ref, ' '.join(v['consensus']))
                else:
                    out += '#\n<msa id="{0}" ref="{1}">\n'.format(k, ref)
                outs = msa2str(v, wordlist=True)
                out += outs
                out += "</msa>\n"

    if distances and 'distances' not in keywords['ignore']:
        out += '\n# DISTANCES\n<dst>\n'
        out += distances + '</dst>\n'

    if trees:
        out += '\n# TREES\n' + trees

    if scorer and 'scorer' not in keywords['ignore']:
        out += '\n# SCORER\n' + scorer

    out += '\n# DATA\n' if keywords['prettify'] else ''
    out += 'ID\t' + '\t'.join(header) + '\n'

    # check for gloss in header to create nice output format
    if formatter in header:
        idx = header.index(formatter)
        formatter = None
        sorted_data = sorted(data.keys(), key=lambda x: data[x][idx])
    elif len(formatter.split(',')) == 2:
        idxA, idxB = formatter.split(',')
        idxA = header.index(idxA)
        idxB = header.index(idxB)
        idx = idxA
        formatter = None
        sorted_data = sorted(data.keys(), key=lambda x: (data[x][idxA], data[x][idxB]))
    else:
        idx = False
        formatter = ''
        sorted_data = sorted(data.keys())

    for key in sorted_data:
        # get the line
        line = data[key]

        # check for formatter
        if idx in range(len(line)):
            if line[idx] != formatter:
                out += '#\n' if keywords['prettify'] else ''
                formatter = line[idx]

        # add the key
        out += text_type(key)

        # add the rest of the values
        for value in line:
            if type(value) == list:
                try:
                    out += '\t' + ' '.join(value)
                except:
                    out += '\t' + ' '.join([text_type(v) for v in value])
            elif type(value) == int:
                out += '\t' + text_type(value)
            elif type(value) == float:
                out += '\t{0:.4f}'.format(value)
            else:
                out += '\t' + value
        out += '\n'

    util.write_text_file(
        filename + '.' + keywords['fileformat'],
        out + keywords.get('stamp', ''),
        normalize="NFC")
    return

def tsv2triple(wordlist, outfile=None):
    """
    Function converts a wordlist to a triple data structure.

    Notes
    -----
    The basic values of which the triples consist are:
      * ID (the ID in the TSV file)
      * COLUMN (the column in the TSV file)
      * VALUE (the entry in the TSV file)
    """
    tstore = []
    for head in wordlist.header:
        log.debug('tsv2triple: ' + head)
        for key in wordlist:
            tstore.append((key, head.upper(), wordlist[key, head]))

    if outfile:
        out = ''
        for a, b, c in tstore:
            if isinstance(c, list):
                c = ' '.join([text_type(x) for x in c])
            if c != '-':
                out += '{0}\t{1}\t{2}\n'.format(a, b, c)
        util.write_text_file(outfile, out, normalize='NFC')
    return tstore


def triple2tsv(triples_or_fname, output="table"):
    """
    Function reads in a triple file and converts it to a tabular data structure.
    """
    D = defaultdict(dict)
    idxs = set()
    cols = set()

    if not isinstance(triples_or_fname, list):
        triples_or_fname = util.read_text_file(
            triples_or_fname, normalize='NFD', lines=True)

    for line in triples_or_fname:
        if isinstance(line, (text_type, str)):
            line = line.split('\t')
        a, b, c = line
        D[a][b] = c
        idxs.add(a)
        cols.add(b)

    idxs = sorted(idxs)
    cols = sorted(cols)
    table = [[idx] + [D.get(idx, {}).get(col, '') for col in cols] for idx in idxs]

    if output not in ['wordlist', 'dict']:
        return [["ID"] + cols] + table

    wlD = {int(line[0]): line[1:] for line in table}
    wlD[0] = cols
    return wlD


def coverage(wordlist):
    """
    Determine the average coverage of a wordlist.
    """
    return {taxon: len(wordlist.get_dict(col=taxon)) for taxon in wordlist.taxa}


def wl2multistate(wordlist, ref):
    """
    Helper function converts a wordlist to multistate format (compatible with PAUP).
    """

    # convert the data to a multistate matrix
    # get etymological dictionary
    wordlist.get_etymdict(ref=ref)

    # define chars, we only have a limited set, unfortunately
    chars = ascii_letters + digits

    # iterate over all cognate sets and assign the chars
    matrix = []
    for c in wordlist.concepts:
        taxon_to_cognate_set = wordlist.get_dict(concept=c, entry='cogid')

        distinct_states = set()
        for taxon in wordlist.taxa:
            distinct_states.update(taxon_to_cognate_set.get(taxon, [0]))

        # make converter
        if len(distinct_states) > len(chars):  # pragma: no cover
            # FIXME: This shouldn't just be a warning, because we will get a KeyError
            # down below, since zip just returns a list of length len(chars)!
            log.warn('more distinct states than available characters!')
        char_map = dict(zip(sorted(distinct_states), chars))
        char_map['-'] = '-'

        line = []
        for taxon in wordlist.taxa:
            states = set(taxon_to_cognate_set.get(taxon, ['-']))
            assert states  # exclude the case len(taxon_to_cognate_set[taxon]) == 0
            if len(states) == 1:
                line.append(char_map[states.pop()])
            else:
                line.append('({0})'.format(
                    "".join([char_map[x] for x in sorted(states)])))

        matrix.append(line)

    return misc.transpose(matrix)
