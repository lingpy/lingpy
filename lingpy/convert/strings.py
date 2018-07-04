# *-* coding: utf-8 *-*
"""
Basic functions for the conversion of Python-internal data into strings.
"""
from __future__ import unicode_literals
import unicodedata
from collections import defaultdict

from lingpy import util
from lingpy.convert.html import template_path

def scorer2str(scorer):
    """
    Convert a scoring function to a string.
    """

    # get sorted representation of characters
    chars = sorted(scorer.chars2int)
    out = ''

    # write stuff to string
    for i, charA in enumerate(chars):
        out += charA
        for j, charB in enumerate(chars):
            out += '\t{0:.2f}'.format(scorer[charA, charB])
        out += '\n'

    return out


def msa2str(msa, wordlist=False, comment="#",
            _arange='{stamp}{comment}\n{meta}{comment}\n{body}', merge=False):
    """
    Function converts an MSA object into a string.
    """

    if 'stamp' in msa:
        stamp = msa['stamp']
    elif hasattr(msa, 'stamp'):
        stamp = msa.stamp
    else:
        stamp = ''

    def mergeit(m, alm, gap='-'):
        """
        Merge a given string according to the merging dictionary m.
        """
        new_alm = ['' for i in range(max(m.values()) + 1)]
        for i, s in enumerate(alm):

            new_alm[m[i]] += s
            if '-' in new_alm[m[i]] and len(new_alm[m[i]]) > 1:
                new_alm[m[i]] = new_alm[m[i]].replace('-', '')

        return new_alm

    meta = ''

    # define global vars for alignment and taxa for convenience
    if wordlist or isinstance(msa, dict):
        taxa = msa['taxa']
        alms = msa['alignment']
    else:
        taxa = msa.taxa
        alms = msa.alignment

    alm_length = len(alms[0])

    # add merge in output as feature
    if hasattr(msa, 'merge') or 'merge' in msa:
        if merge:
            try:
                merger = msa['merge']
            except:
                merger = msa.merge

            new_alms = []

            for alm in alms:
                new_alms += [mergeit(merger, alm)]
            alms = new_alms
        else:
            merger = dict([(i, i) for i in range(len(alms[0]))])
    else:
        merger = dict([(i, i) for i in range(len(alms[0]))])

    body = ''

    # if wordlist ist set to True, don't write the header line and put the
    # after comment
    if wordlist:
        # get formatter
        formatter = max([len(t) for t in msa['taxa'] + ['COLUMNID']])
        for a, b, c in zip(msa['ID'], taxa, alms):
            body += '{0}\t{1}'.format(a, b.ljust(formatter, '.')) + '\t'
            body += '\t'.join(c) + '\n'
        alm_len = len(c)

    elif isinstance(msa, dict):
        # get formatter
        formatter = max([len(t) for t in msa['taxa']])
        body += msa['dataset'] + '\n'
        body += msa['seq_id'] + '\n'
        for a, b in zip(taxa, alms):
            body += a.ljust(formatter, '.') + '\t'
            body += '\t'.join(b) + '\n'
        alm_len = len(b)
    else:
        # get formatter
        formatter = max([len(t) for t in msa.taxa])
        body += msa.dataset + '\n'
        body += msa.seq_id + '\n'
        for a, b in zip(taxa, alms):
            body += a.ljust(formatter, '.') + '\t'
            body += '\t'.join(b) + '\n'
        alm_len = len(b)

    if 'local' in msa:
        local = msa['local']
    elif hasattr(msa, 'local'):
        local = msa.local
    else:
        local = False

    if 'swaps' in msa:
        swaps = msa['swaps']
    elif hasattr(msa, 'swaps'):
        swaps = msa.swaps
    else:
        swaps = False

    if 'consensus' in msa:
        consensus = msa['consensus']
    elif hasattr(msa, 'consensus'):
        consensus = msa.consensus
    else:
        consensus = False

    if wordlist:
        meta += '0\t' + 'COLUMNID'.ljust(formatter, '.') + '\t' + '\t'.join([
            str(i + 1) for i in range(alm_len)
        ])
        meta += '\n#\n'

    if local:
        if wordlist:
            meta += '{0}\t{1}\t'.format(0, 'LOCAL'.ljust(formatter, '.'))
        else:
            meta += '{0}\t'.format('LOCAL'.ljust(formatter, '.'))
        tmp = []
        for i in range(alm_len):
            if i in local:
                tmp += ['*']
            else:
                tmp += ['.']
        meta += '\t'.join(mergeit(merger, tmp)) + '\n'
    if swaps:
        if wordlist:
            meta += '{0}\t{1}\t'.format(0, 'CROSSED'.ljust(formatter, '.'))
        else:
            meta += '{0}\t'.format('SWAPS'.ljust(formatter, '.'))
        tmp = alm_length * ['.']
        for swap in swaps:
            a, b, c = swap
            tmp[a] = '+'
            tmp[b] = '-'
            tmp[c] = '+'
        meta += '\t'.join(mergeit(merger, tmp, gap='')) + '\n'

    if consensus:
        if wordlist:
            meta += '{0}\t{1}\t'.format(0, 'CONSENSUS'.ljust(formatter, '.'))
        else:
            meta += '{0}\t'.format('CONSE'.ljust(formatter, '.'))
        meta += '\t'.join(mergeit(merger, consensus, gap='')) + '\n'

    return _arange.format(
        stamp=stamp,
        meta=meta,
        body=body,
        comment=comment
    )


def matrix2dst(
    matrix,
    taxa=None,
    stamp='',
    filename='',
    taxlen=10,
    comment='#'
):
    """
    Convert matrix to dst-format.

    Parameters
    ----------
    taxa : {None, list}
        List of taxon names corresponding to the distances. Make sure that you
        only use alphanumeric characters and the understroke for assigning the
        taxon names. Especially avoid the usage of brackets, since this will
        confuse many phylogenetic programs.
    stamp : str (default='')
        Convenience stamp passed as a comment that can be used to indicate how
        the matrix was created.
    filename : str
        If you specify a filename, the data will be written to file.
    taxlen : int (default=10)
        Indicate how long the taxon names are allowed to be. The Phylip package
        only allows taxon names consisting of maximally 10 characters. Other
        packages, however, allow more. If Phylip compatibility is not important
        for you and you just want to allow for as long taxon names as possible,
        set this value to 0.
    comment : str (default = '#')
        The comment character to be used when adding additional information in
        the "stamp".

    Returns
    -------
    output : {str or file}
        Depending on your settings, this function returns a string in DST
        (=Phylip) format, or a file containing the string.

    """
    if not taxa:
        taxa = ['t_{0}'.format(i + 1) for i in range(len(matrix))]

    out = ' {0}\n'.format(len(taxa))
    for i, taxon in enumerate(taxa):

        # check for zero-taxlen
        if taxlen == 0:
            dummy = '{0}\t'
            idx = len(taxon)
            joinchar = '\t'  # normally in Phylip this is a space
        else:
            dummy = '{0:' + str(taxlen) + '}'
            idx = taxlen + 1
            joinchar = ' '

        out += dummy.format(taxon)[:idx] + joinchar
        out += joinchar.join(['{0:.2f}'.format(d) for d in
                              matrix[i]])
        out += '\n'
    if stamp:
        out += '{1} {0}'.format(stamp, comment)
    if not filename:
        return out
    else:
        util.write_text_file(filename + '.dst', out)


def pap2nex(
    taxa,
    paps,
    missing=0,
    filename='',
    datatype='STANDARD'
):
    """
    Function converts a list of paps into nexus file format.

    Parameters
    ----------
    taxa : list
        List of taxa.
    paps : {list, dict}
        A two-dimensional list with the first dimension being identical to the
        number of taxa and the second dimension being identical to the number
        of paps. If a dictionary is passed, each key represents a given pap.
        The following two structures will thus be treated identically::

          >>> paps = [[1,0],[1,0],[1,0]] # two languages, three paps
          >>> paps = {1:[1,0], 2:[1,0], 3:[1,0]} # two languages, three paps

    missing : {str, int} (default=0)
        Indicate how missing characters are represented in the original data.

    """
    out = '#NEXUS\n\nBEGIN DATA;\nDIMENSIONS ntax={0} NCHAR={1};\n'
    out += "FORMAT DATATYPE={5} GAP=- MISSING={2} interleave=yes;\n"
    out += "MATRIX\n\n{3}\n;\n\nEND;\n"
    out += "[PAPS-REFERENCE]\n{4}"

    # get longest taxon
    maxTax = max([len(taxon) for taxon in taxa])

    # check whether paps are dict or list
    if hasattr(paps, 'keys'):
        new_paps = [paps[k] for k in sorted(paps)]
        reference = [k for k in sorted(paps)]
    else:
        new_paps = paps
        reference = [k for k in range(1, len(paps) + 1)]

    # create reference
    ref_string = ''
    for i, ref in enumerate(reference):
        ref_string += '[{0} :: {1}]\n'.format(i, ref)
    # create the matrix
    matrix = ""

    for i, taxon in enumerate(taxa):
        tmp = '{0:XXX} '
        matrix += tmp.replace('XXX', str(maxTax)).format(taxon)
        matrix += ''.join([str(itm[i]) for itm in new_paps])
        matrix += '\n'

    if not filename:
        return out.format(
            len(taxa),
            len(paps),
            missing,
            matrix,
            ref_string,
            datatype
        )
    util.write_text_file(
        filename + '.nex',
        out.format(len(taxa), len(paps), missing, matrix, ref_string, datatype))
    return


def pap2csv(
    taxa,
    paps,
    filename=''
):
    """
    Write paps created by the Wordlist class to a csv-file.
    """

    out = "ID\t" + '\t'.join(taxa) + '\n'

    for key in sorted(paps):
        out += '{0}\t{1}\n'.format(
            key,
            '\t'.join(str(i) for i in paps[key])
        )

    if not filename:
        return out
    util.write_text_file(filename + '.csv', out)
    return


def multistate2nex(taxa, matrix, filename='', missing="?"):
    """
    Convert the data in a given wordlist to NEXUS-format for multistate analyses in PAUP.

    Parameters
    ----------
    taxa : list
        The list of taxa that shall be written to file.
    matrix : list
        The multi-state matrix with the first dimension indicating the taxa,
        and the second their states.
    filename : str (default="")
        If not specified, the filename of the Wordlist will be taken,
        otherwise, it specifies the name of the file to which the data will be
        written.
    """

    # set up the nexus template
    nexus = """#NEXUS

BEGIN DATA;
DIMENSIONS ntax={ntax} NCHAR={nchar};
FORMAT RESPECTCASE DATATYPE=STANDARD symbols="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOP0123456789" GAP=? MISSING={missing} interleave=yes;
OPTIONS MSTAXA = POLYMORPH;

MATRIX

{matrix}

END;
"""

    # calculate maximal length of taxon strings
    tlen = max([len(t) for t in taxa])

    # calculate the matrix-text in the nexus template
    matrix_text = ""
    for taxon, line in zip(taxa, matrix):
        ntaxon = taxon + tlen * ' ' + ' '
        ntaxon = ntaxon[:tlen]
        matrix_text += "{0} {1}\n".format(ntaxon, ''.join(line))

    if filename:
        util.write_text_file(
            filename,
            nexus.format(
                ntax=len(taxa),
                nchar=len(matrix[0]),
                matrix=matrix_text,
                missing=missing
            )
        )
    else:
        raise ValueError("[!] A wrong filename was specified!")
    return


def _is_constant(e):
    """
    Returns whether an etymdict dictionary is constant. Helper function for
    `write_nexus`
    
    Parameters
    ----------
    e : dict
        An etymdict sub-dictionary.
    
    Returns
    -------
    bool : bool
        A boolean True/False for whether the character is constant or not.
    """
    if all([_ != 0 for _ in e]):  # constant present
        return True
    elif all([_ == 0 for _ in e]):  # constant absent
        return True
    return False


def write_nexus(
        wordlist,
        mode='mrbayes',
        filename="mrbayes.nex",
        ref="cogid",
        missing="?", gap="-",
        custom=None,
        custom_name='lingpy',
        commands=None, commands_name="mrbayes"):
    """Write a nexus file for phylogenetic analyses.

    Parameters
    ----------
    wordlist : lingpy.basic.wordlist.Wordlist
        A Wordlist object containing cognate IDs.
    mode : str (default="mrbayes")
        The name of the output nexus style. Valid values are:
            * 'MRBAYES': a MrBayes formatted nexus file.
            * 'SPLITSTREE': a SPLITSTREE formatted nexus file.
            * 'BEAST': a BEAST formatted nexus file.
            * 'BEASTWORDS': a BEAST formatted nexus for word-partitioned
               analyses.
            * 'TRAITLAB': a TRAITLab formatted nexus.
    filename : str (default=None)
        Name of the file to which the nexus file will be written.
        If set to c{None}, then this function will not write the nexus ontent
        to a file, but simply return the content as a string.
    ref: str (default="cogid")
        Column in which you store the cognate sets in your data.
    gap : str (default="-")
        The symbol for gaps (not relevant for linguistic analyses).
    missing : str (default="?")
        The symbol for missing characters.
    custom : list {default=None)
        This information allows to add custom information to the nexus file, like, for
        example, the structure of the characters, their original concept, or their
        type, and it will be written into a custom block in the nexus file. The name of
        the custom block can be specified with help of the `custom_name` keyword. The
        content is a list of strings which will be written line by line into the custom
        block.
    custom_name : str (default="lingpy")
        The name of the custom block which will be written to the file.
    commands : list (default=None)
        If specified, will write an additional block containing commands for
        phylogenetic software. The commands are passed as a list, containing
        strings. The name of the block is given by the keywords commands_name.
    commands_name : str (default="mrbayes")
        Determines how the block will be called to which the commands will be
        written.

    Returns
    -------
    nexus : str
        A string containing nexus file output
    """
    templates = {
        'BEAST': 'beast.nex',
        'BEASTWORDS': 'beast.nex',
        'SPLITSTREE': 'splitstree.nex',
        'MRBAYES': 'mrbayes.nex',
        'TRAITLAB': 'splitstree.nex',
    }
    
    block = "\n\nBEGIN {0};\n{1}\nEND;\n"  # template for nexus blocks
    
    # check for valid mode
    mode = mode.upper()
    if mode not in templates.keys():
        raise ValueError("Unknown output mode %s" % mode)

    # check for valid template
    template = templates.get(mode)
    tpath = util.Path(template_path(template))
    if tpath.exists:
        _template = util.read_text_file(tpath.as_posix())
    else:  # pragma: no cover
        raise IOError("Unknown template %s" % template)

    # check that `ref` is a valid column
    if ref not in wordlist._alias:
        raise KeyError("Unknown _ref_ column in wordlist '%s'" % ref)

    # retrieve the matrix
    matrix = [[] for x in range(wordlist.width)]
    etd = wordlist.get_etymdict(ref=ref)
    concepts = sorted([(cogid, wordlist[[
        x[0] for x in vals if x][0]][wordlist._rowIdx]) for (cogid, vals) in
        etd.items()],
        key=lambda x: (x[1], x[0]))
    # and missing data..
    missing_ = {t: [concept for (cogid, concept) in concepts if concept not in wordlist.get_list(
                col=t, entry=wordlist._row_name, flat=True)] for t in
                wordlist.cols}
    
    # add ascertainment character for mode=BEAST
    if mode == 'BEAST':
        matrix = [['0'] for m in matrix]
    
    # skip the constant sites for traitlab
    if mode == 'TRAITLAB':
        concepts = [(i, c) for (i, c) in concepts if not _is_constant(etd[i])]
    
    # fill matrix
    for i, t in enumerate(wordlist.cols):
        previous = ''
        for cogid, concept in concepts:
            if previous != concept:
                previous = concept
                # add ascertainment character for mode=BEASTWORDS. Note that if
                # a given word:language is missing, then its ascertainment
                # character is the `missing` character.
                if mode == "BEASTWORDS":
                    matrix[i] += ['0'] if concept not in missing_[t] else [missing]
            matrix[i] += ['1'] if etd[cogid][i] else ['0'] if concept not in \
                missing_[t] else [missing]

    # parse characters into `charsets` (a dict of word=>siteindex positions),
    # and `chars` (a list of characters).
    charsets, chars, previous = defaultdict(list), [], ''
    for i, (cogid, concept) in enumerate(concepts, 1):
        char = util.nexus_slug(concept)
        # add label for ascertainment character in BEAST mode
        if i == 1 and mode == 'BEAST':
            chars.append("_ascertainment")
        # add label for per-word ascertainment characters in BEASTWORDS
        if mode == 'BEASTWORDS' and previous != concept:
            chars.append("%s_ascertainment" % char)
            charsets[char].append(len(chars))
        # finally add label.
        chars.append(char)
        charsets[char].append(len(chars))
        previous = concept
    
    # create character labels block if needed
    if mode in ('BEAST', 'BEASTWORDS'):
        charblock = ",\n".join(["\t%d %s" % o for o in enumerate(chars, 1)])
    else:
        charblock = ""
    
    # create charsets block
    blockname, assumptions = None, ""
    if mode in ('BEASTWORDS', 'MRBAYES'):
        charsets = ["\tcharset %s = %d-%d;" % (
            c, min(m), max(m)) for (c, m) in charsets.items()
        ]
        blockname = 'ASSUMPTIONS' if mode == 'BEASTWORDS' else 'MRBAYES'
        assumptions = "\n".join(charsets)
    
    # commands
    if commands_name.upper() == blockname and len(assumptions) and commands:
        # merge commands specified in function call into output blockname
        assumptions += "\n" + "\n".join("\t%s" % c for c in commands)
    else:
        # different commands block set in commands_name.
        assumptions += block.format(commands_name, '\n'.join(commands)) if commands else ''
    
    # convert state matrix to string.
    _matrix = ""
    maxtaxlen = max([len(util.nexus_slug(t)) for t in wordlist.cols]) + 1
    for i, (taxon, m) in enumerate(zip(wordlist.cols, matrix)):
        _matrix += str(util.nexus_slug(taxon) + maxtaxlen * ' ')[:maxtaxlen] + ' '
        _matrix += ''.join([
            '({0})'.format(c) if len(c) > 1 else str(c) for c in m
        ]) + '\n'
    _matrix = _matrix.rstrip()  # remove trailing

    # TODO: symbols could be more than "01" but we this function doesn't handle
    # multistate data so we just specify them here.
    symbols = '01'

    text = _template.format(
        matrix=_matrix,
        ntax=wordlist.width,
        nchar=len(matrix[0]),
        gap=gap, missing=missing,
        dtype='RESTRICTION' if mode == 'MRBAYES' else 'STANDARD',
        commands=block.format(blockname, assumptions),
        custom=block.format(custom_name, '\n'.join(custom)) if custom else '',
        symbols=symbols, chars=charblock
    )
    text = text.replace("\t", " " * 4)  # normalise tab-stops
    for i, (cogid, concept) in enumerate(concepts, 1):
        text += '\n[MATRIX:{0}=COGID:{1}=CONCEPT:{2}]'.format(i, cogid, concept)
    if filename:
        util.write_text_file(filename, text)
    return text
