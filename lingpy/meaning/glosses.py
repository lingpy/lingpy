"""
Module provides functions for the handling of concept glosses in linguistic datasets.
"""
import re

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
            with each tuple consisting of 5 items, namely:
        * the main part ("main"),
        * the start character indicating a potential comment ("comment_start"),
        * the comment (everything occurring in brackets in the input string
          ("comment"),
        * the end character indicating the end of a potential comment
          ("comment_end"),
        * the part of speech, in case this was specificied by a preceding
          "the" or a preceding "to" in the mainpart of the string ("pos"),
        * the prefix, that is, words, like, eg. "be", "in", which may precede the main
          gloss in concept lists, as in "be quiet" ("prefix")

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
    constituents = re.split(',|;', gloss)

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

            G += [(mainpart, comment_start, comment, comment_end, pos, prefix)]
    
    if output == 'dict':
        return [dict(zip(["main","comment_start", "comment", "comment_end",
            "pos", "prefix"], g)) for g in G]
    return G

if __name__ == '__main__':
    for word in ['to kill (verb)', 'kill [somebody]', 'work, the ']:
        print(parse_gloss(word))

