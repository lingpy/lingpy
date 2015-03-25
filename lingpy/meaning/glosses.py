"""
Module provides functions for the handling of concept glosses in linguistic datasets.
"""
import re
from ..read.csv import csv2list
from ..align.pairwise import edit_dist

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
        A list of tuples, or a list of dictionaries, 
            with each tuple consisting of 6 items, namely:
        * the main part ("main"),
        * the start character indicating a potential comment ("comment_start"),
        * the comment (everything occurring in brackets in the input string
          ("comment"),
        * the end character indicating the end of a potential comment
          ("comment_end"),
        * the part of speech, in case this was specificied by a preceding
          "the" or a preceding "to" in the mainpart of the string ("pos"),
        * the prefix, that is, words, like, eg. "be", "in", which may precede the main
          gloss in concept lists, as in "be quiet" ("prefix"),
        * the longest constituent, which is identical with the main part if
            there's no whitespace in the main part, otherwise the longest part
            part of the main gloss split by whitespace ("longest_part")
        * the original gloss (for the purpose of testing, labelled "gloss")

        If "dict" is chosen as output, this returns a list of dictionaries with
        the keys as specified in brackets above.

    Note
    ----
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
            mainpart = ''.join([m for m in mainpart if m not in '?!"¨:;,»«´“”*'])
            mainpart = mainpart.strip()
            mainpart = mainpart.lower()

            # search for pos-markers
            pos_markers = [
                    ('the', 'noun'),
                    ('a','noun'),
                    ('to','verb'),
                    ]
            pos = ''

            for p,t in pos_markers:
                if mainpart.startswith(p+' ') and not pos:
                    mainpart = mainpart[len(p)+1:]
                    pos = t
            
            # search for strip-off-prefixes
            prefixes = ['be', 'in', 'at']
            prefix = ''
            for p in prefixes:
                if mainpart.startswith(p+' ') and not prefix:
                    mainpart = mainpart[len(p)+1:]
                    prefix = p
            
            # search for pos in comment
            abbreviations = [
                    ('vb','verb'),
                    ('v.','verb'),
                    ('v','verb'),
                    ('adj','adjective'),
                    ('nn','noun'),
                    ('n.','noun'),
                    ('adv','adverb')
                    ]
            for p,t in sorted(abbreviations, key=lambda x: len(x),
                    reverse=True):
                if p in comment and not pos:
                    pos = t             

            # check for a "first part" in case we encounter white space in the
            # data (and return only the largest string of them)
            if ' ' in mainpart:
                parts = mainpart.split(' ')
                best_part = sorted(parts, key=lambda x: len(x))[-1]
            else:
                best_part = mainpart

            G += [(mainpart, comment_start, comment, comment_end, pos, prefix,
                best_part, constituent)]
    
    if output == 'dict':
        return [dict(zip(["main","comment_start", "comment", "comment_end",
            "pos", "prefix", "longest_part", "gloss"], g)) for g in G]
    return G


def compare_conceptlists(list1, list2, output='', match=None):
    """
    Function compares two concept lists and outputs suggestions for mapping.

    Notes
    -----
    Idea is to take one conceptlist as the basic list and then to search for a
    plausible mapping of concepts in the second list to the first list. All
    suggestions can then be output in various forms, both with multiple matches
    excluded or included, and in textform or in other forms.
    """
    
    # check for match quality
    if not match:
        match = [1,2,3,4]
        
    # take first list as basic list
    base = csv2list(list1)
    comp = csv2list(list2)

    # get headers
    baseh,base = base[0],base[1:]
    comph,comp = comp[0],comp[1:]
    
    # make sure to raise if 'gloss' is not in the headers
    if (not "GLOSS" in baseh and not "GLOSS" in comph) or \
            (not "NUMBER" in baseh and not "NUMBER" in comph):
        raise ValueError("There is not field for 'GLOSS' or 'NUMBER' in the header of the input lists.")
    
    # get gloss indices
    bidx = baseh.index('GLOSS')
    cidx = comph.index('GLOSS')
    bnum = baseh.index("NUMBER")
    cnum = comph.index("NUMBER")
    
    # extract glossing information from the data
    B = {}
    idx = 1
    for i,line in enumerate(base):
        gloss = line[bidx]
        gdata = parse_gloss(gloss, output='dict')
        for gdatum in gdata:
            gdatum['number'] = line[bnum] # we won't need "enumerate" XXX
            B[idx] = gdatum
            idx += 1
    
    idx = 1
    C = {}
    for i,line in enumerate(comp):
        gloss = line[cidx]
        gdata = parse_gloss(gloss, output = 'dict')
        for gdatum in gdata:
            gdatum['number'] = line[cnum] # we won't need "enumerate" XXX
            C[idx] = gdatum
            idx += 1

    # now that we have prepared all the glossed list as planned, we compare
    # them item by item and check for similarity
    sims = []
    for i,a in sorted(B.items()):
        for j,b in sorted(C.items()):

            # first-order-match: identical glosses
            if a['gloss'] == b['gloss']:
                sims += [(i,j,1)]
            # second-order match: identical main-parts
            elif a['main'] == b['gloss'] or a['gloss'] == b['main'] or \
                    a['main'] == b['main']:
                # best match if pos matches 
                if a['pos'] == b['pos']:
                    sims += [(i,j,2)]

                # less good match if pos mismatches
                else:
                    sims += [(i,j,3)]
            elif a['longest_part'] == b['longest_part']:
                sims += [(i,j,4)]
            elif b['longest_part'] in a['gloss']:
                sims += [(i,j,5)]
            elif a['longest_part'] in b['gloss']:
                sims += [(i,j,6)]
    
    # get the number of items which were not matched in the second list
    matched = [x[1] for x in sims if x[2] in match]
    not_matched = [idx for idx in C if idx not in matched]
    for idx in not_matched:
        sims += [(0,idx,7)]

    # prepare the output
    output = []
    for a,b,c in sims:
        if c in match:
            output += [(c,B[a]['gloss'],B[a]['number'], C[b]['gloss'],C[b]['number'])]
        elif c not in match and c < 7:
            pass
        else:
            output += [(c,'?','?',C[b]['gloss'],C[b]['number'])]
            
    return output

    
