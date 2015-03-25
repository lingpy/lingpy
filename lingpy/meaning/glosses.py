"""
Module provides functions for the handling of concept glosses in linguistic datasets.
"""
import re

def parse_gloss(gloss):
    """
    Parse a gloss into its constituents by applying some general logic.

    Parameters
    ----------
    gloss : str
        The gloss as found in various sources (we assume that we are dealing
        with English glosses here.

    Returns
    -------
    constituents : list
        A list of tuples, with each tuple consisting of 4 items, namely:
        * the main part,
        * the start character indicating a potential comment,
        * the comment (everything occurring in brackets in the input string,
        * the end character indicating the end of a potential comment, and
        * the part of speech, in case this was specificied by a preceding
          "the" or a preceding "to" in the mainpart of the string.

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
    
    return G

if __name__ == '__main__':
    for word in ['to kill (verb)', 'kill [somebody]', 'work, the ']:
        print(parse_gloss(word))

