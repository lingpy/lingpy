# *-* coding: utf-8 *-*
"""
Module provides functions for the handling of concept glosses in linguistic datasets.
"""
from __future__ import print_function, division, unicode_literals
import re

from lingpy.read.csv import csv2list
from lingpy.util import write_text_file


def parse_gloss(gloss, output='list'):
    """
    Parse a gloss into its constituents by applying some general logic.

    Parameters
    ----------
    gloss : str
        The gloss as found in various sources (we assume that we are dealing
        with English glosses here.
    output : str (default="list")
        Determine the output of the parsing routine. Select between "list"
        which will return a list of tuples, or "dict", which will return a list
        of dictionaries.

    Returns
    -------
    constituents : {list, dict}
        A list of tuples, or a list of dictionaries, with each tuple consisting of 6 items, namely:
        * the main part ("main"),
        * the start character indicating a potential comment ("comment_start"),
        * the comment (everything occurring in brackets in the input string ("comment"),
        * the end character indicating the end of a potential comment ("comment_end"),
        * the part of speech, in case this was specificied by a preceding "the" or a preceding "to" in the mainpart of the string ("pos"),
        * the prefix, that is, words, like, eg. "be", "in", which may precede the main gloss in concept lists, as in "be quiet" ("prefix"),
        * the longest constituent, which is identical with the main part if there's no whitespace in the main part, otherwise the longest part part of the main gloss split by whitespace ("longest_part")
        * the parts of a gloss, if the constituent contains multiple words ("parts")
        * the original gloss (for the purpose of testing, labelled "gloss")
        
        If "dict" is chosen as output, this returns a list of dictionaries with
        the keys as specified in brackets above.

    Notes
    -----

    The basic purpose of this function is to provide a means to make it easier
    to compare meanings across different resources. Often, linguists will
    annotate their resources quite differently, and for one and the same
    concept, we may find very different glosses. The concept "kill [verb]", for
    example may be glossed as "to kill", "kill", "kill (v.)", "kill
    (somebody)", etc. In order to guarantee comparability, this function tries
    to use basic knowledge of glossing tendencies to disentangle the variety of
    glossing styles which can be found in the literature. Thus, in the case of
    "kill [verb]", the function will analyze the different strings as follows::

        >>> from lingpy.meaning.glosses import parse_gloss
        >>> glosses = ["to kill", "kill", "kill (v.)", "kill (somebody)"]
        >>> for gloss in glosses:
        ...     parsed_gloss = parse_gloss(gloss, output='dict')
        ...     print(parsed_gloss[0]['gloss'], parsed_gloss[0]['pos'])
        kill verb
        kill
        kill verb
        kill

    As can be seen: it seeks to extract the most important part of the gloss
    and may thus help to compare different glosses across different resources.
    """
    # check for constituents
    constituents = re.split(',|;|/| or ', gloss)

    # create output dictionary
    G = []

    gpos = ''

    for constituent in constituents:
        if constituent.strip():
            # separate brackets as comments
            mainpart = ''
            comment = ''
            comment_start = ''
            comment_end = ''
            commented = False
            for char in constituent:

                if commented:
                    comment += char
                else:
                    mainpart += char

                # check for commented regions, we use some basic separators here
                if char in '([{（<':
                    commented = True
                    comment_start += char
                    mainpart = mainpart[:-1]
                if char in ')]}）>':
                    commented = False
                    comment_end += char
                    comment = comment[:-1]

            # search for useless chars
            mainpart = ''.join([m for m in mainpart if m not in '?!"¨:;,»«´“”*+-'])
            mainpart = mainpart.strip()
            mainpart = mainpart.lower()

            # search for pos-markers
            pos_markers = [('the', 'noun'), ('a', 'noun'), ('to', 'verb')]
            pos = ''

            if gpos:
                pos = gpos
            else:
                for p, t in pos_markers:
                    if mainpart.startswith(p + ' ') and not pos:
                        mainpart = mainpart[len(p) + 1:]
                        pos = t
                        gpos = t

            # search for strip-off-prefixes
            prefixes = ['be', 'in', 'at']
            prefix = ''
            for p in prefixes:
                if mainpart.startswith(p + ' ') and not prefix:
                    mainpart = mainpart[len(p) + 1:]
                    prefix = p

            # check for a "first part" in case we encounter white space in the
            # data (and return only the largest string of them)
            if ' ' in mainpart:
                parts = mainpart.split(' ')
                best_part = sorted(parts, key=lambda x: len(x))[-1]
            else:
                best_part = mainpart
                parts = [mainpart]

            # search for pos in comment
            abbreviations = [
                ('vb', 'verb'),
                ('v.', 'verb'),
                ('v', 'verb'),
                ('adj', 'adjective'),
                ('nn', 'noun'),
                ('n.', 'noun'),
                ('adv', 'adverb'),
                ('noun', 'nount'),
                ('verb', 'verb'),
                ('adjective', 'adjective'),
                ('cls', 'classifier')
            ]
            for p, t in sorted(abbreviations, key=lambda x: len(x), reverse=True):
                if p in comment and not pos:
                    pos = t
                elif p in parts:
                    pos = t
                elif t in comment and not pos:
                    pos = t
                elif t in parts:
                    pos = t

            G += [(
                mainpart,
                comment_start,
                comment,
                comment_end,
                pos,
                prefix,
                best_part,
                parts,
                gloss)]

    if output == 'dict':
        return [dict(zip([
            "main",
            "comment_start",
            "comment",
            "comment_end",
            "pos",
            "prefix",
            "longest_part",
            "parts",
            "gloss"
        ], g)) for g in G]
    return G


def compare_concepts(c1, c2):
    """
    Debug-function for concept comparison.
    """
    c1g = parse_gloss(c1, output='dict')
    c2g = parse_gloss(c2, output='dict')

    sims = []
    for a in c1g:
        for b in c2g:
            # first-order-match: identical glosses
            if a['gloss'] == b['gloss']:
                sims += [(a['gloss'], b['gloss'], 1)]  # [(i,j,1)]
            # second-order match: identical main-parts
            elif a['main'] == b['gloss'] or a['gloss'] == b['main'] or \
                    a['main'] == b['main']:
                # best match if pos matches
                if a['pos'] == b['pos'] and a['pos']:
                    sims += [(a['gloss'], b['gloss'], 2)]  # [(i,j,2)]
                # less good match if pos mismatches
                else:
                    sims += [(a['gloss'], b['gloss'], 3)]  # [(i,j,3)]
            elif a['longest_part'] == b['longest_part']:
                if a['pos'] == b['pos'] and a['pos']:
                    sims += [(a['gloss'], b['gloss'], 4)]  # [(i,j,4)]
                else:
                    sims += [(a['gloss'], b['gloss'], 5)]  # [(i,j,5)]
            elif b['longest_part'] in a['parts']:
                sims += [(a['gloss'], b['gloss'], 6)]  # [(i,j,6)]
            elif a['longest_part'] in b['parts']:
                sims += [(a['gloss'], b['gloss'], 7)]  # [(i,j,7)]
            elif not sims:
                sims += [('?', '?', 8)]
    return sims

def compare_conceptlists(
        list1,
        list2,
        output='',
        match=None,
        filename='matches',
        **keywords):
    """
    Function compares two concept lists and outputs suggestions for mapping.

    Notes
    -----
    Idea is to take one conceptlist as the basic list and then to search for a
    plausible mapping of concepts in the second list to the first list. All
    suggestions can then be output in various forms, both with multiple matches
    excluded or included, and in textform or in other forms.

    What is important, regarding the output here, is, that the output contains
    all matches, including non-matched items which occur **in the second list
    but not in the first list**. Non-matched items which occur in the first
    list but not in the second list are ignored.

    The syntax for matching types is organized as follows:

    * 1 indicates a full match between glosses, including information on part
      speech and the like
    * 2 indicates a very good match between a full gloss and the main part of a
      gloss or the two main parts of a gloss
    * 3 indicates a very good match between the main parts of two glosses with
      non-matching information regarding part of speech
    * 4 indicates that the longest part of two glosses matches along with the
      part-of-speech information.
    * 5 indicates that the longest part of two glosses matches with
      non-matching part-of-speech information.
    * 6 indicates that the longest part of the first list is matched by one of
      the parts in the second list
    * 7 indicates that the longest part of the second list is matched by one of
      the parts in the first list
    * 8 indicates that no match could be found.
    """
    # check for match quality
    if not match:
        match = [1, 2, 3, 4, 5]

    # check for keywords
    defaults = dict(
        id_name='CONCEPTICON_ID',
        gloss_name='CONCEPTICON_GLOSS',
        match_quality='MATCH_QUALITY',
        gloss='GLOSS',
        number='NUMBER')
    defaults.update(keywords)

    # take first list as basic list
    base = csv2list(list1)
    comp = csv2list(list2)

    # get headers
    baseh, base = base[0], base[1:]
    comph, comp = comp[0], comp[1:]

    # make sure to raise if 'gloss' is not in the headers
    if (not defaults["gloss"] in baseh and not defaults["gloss"] in comph) or \
            (not defaults["number"] in baseh and not defaults["number"] in comph):
        raise ValueError(
            "[!] There is no field for '{0}' or '{1}'".format(
                keywords['gloss'],
                keywords['number']
            ) + " in the header of the input lists.")

    # get gloss indices
    bidx = baseh.index(defaults['gloss'])
    cidx = comph.index(defaults['gloss'])
    bnum = baseh.index(defaults['number'])
    cnum = comph.index(defaults['number'])

    # extract glossing information from the data
    B = {}
    idx = 1
    for i, line in enumerate(base):
        gloss = line[bidx]
        gdata = parse_gloss(gloss, output='dict')
        for gdatum in gdata:
            gdatum['number'] = line[bnum]  # we won't need "enumerate" XXX
            B[idx] = gdatum
            idx += 1

    idx = 1
    line2idx = {}
    C = {}
    for i, line in enumerate(comp):
        gloss = line[cidx]
        gdata = parse_gloss(gloss, output='dict')
        for gdatum in gdata:
            gdatum['number'] = line[cnum]  # we won't need "enumerate" XXX
            C[idx] = gdatum
            try:
                line2idx[i] += [idx]
            except KeyError:
                line2idx[i] = [idx]
            idx += 1

    # now that we have prepared all the glossed list as planned, we compare
    # them item by item and check for similarity
    sims = []
    for i, a in sorted(B.items()):
        for j, b in sorted(C.items()):
            # first-order-match: identical glosses
            if a['gloss'] == b['gloss']:
                sims += [(i, j, 1)]
            # second-order match: identical main-parts
            elif a['main'] == b['gloss'] or a['gloss'] == b['main'] or \
                    a['main'] == b['main']:
                # best match if pos matches
                if a['pos'] == b['pos']:
                    sims += [(i, j, 2)]
                # less good match if pos mismatches
                else:
                    sims += [(i, j, 3)]
            elif a['longest_part'] == b['longest_part']:
                if a['pos'] == b['pos'] and a['pos']:
                    sims += [(i, j, 4)]
                else:
                    sims += [(i, j, 5)]
            elif b['longest_part'] in a['parts']:
                sims += [(i, j, 6)]
            elif a['longest_part'] in b['parts']:
                sims += [(i, j, 7)]

    # get the number of items which were not matched in the second list
    matched = [x[1] for x in sims if x[2] in match]
    not_matched = [idx_ for idx_ in C if idx_ not in matched]
    for idx in not_matched:
        sims += [(0, idx, 8)]

    # sort the matches, add them to a dictionary
    best = {}
    for a, b, c in sims:
        try:
            best[b] += [(a, c)]
        except KeyError:
            best[b] = [(a, c)]

    for k, v in best.items():
        best[k] = sorted(set(v), key=lambda x: x[1])

        if best[k][0][1] in matched:
            best[k] = [best[k][0]]

    # prepare the output
    out = []
    for b in best:  # in sims:
        for a, c in best[b]:
            if c in match:
                out += [(c, B[a]['gloss'], B[a]['number'], C[b]['gloss'], C[b]['number'])]
            elif c == 0:
                out += [(c, '?', '0', C[b]['gloss'], C[b]['number'])]

    if not output:
        return out
    elif output == 'tsv':
        added = []
        txt = ['\t'.join(comph) + '\t{0}\t{1}\t{2}\n'.format(
            defaults['id_name'],
            defaults['gloss_name'],
            defaults['match_quality'])]
        for i, line in enumerate(comp):
            for idx in line2idx[i]:
                if idx in best:
                    data = best[idx]
                else:
                    data = [('?', '0')]

                for a, b in data:
                    if b in match or b == 8:
                        try:
                            base_gloss = B[a]['gloss']
                            base_num = B[a]['number']
                        except KeyError:
                            base_gloss = '???'
                            base_num = '0'

                        nline = '\t'.join(line) + '\t' + str(base_num) + '\t' + \
                                base_gloss + '\t' + str(b) + '\n'
                        if nline not in added:
                            txt += [nline]
                            added += [nline]
                    else:
                        nline = '\t'.join(line) + '\t???\t???\t8\n'
                        if nline not in added:
                            txt += [nline]
                            added += [nline]

            txt[-1] += '\n'

        out = [txt[0]] + sorted(txt[1:], key=lambda x: x[x.index('\t')])
        write_text_file(filename, ''.join(out))
