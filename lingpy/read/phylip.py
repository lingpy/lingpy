"""
Module provides functions to read in various formats from the Phylip package.
"""
import regex as re

def read_dst(
        filename
        ):
    """
    Function reads files in Phylip dst-format.

    Parameters
    ----------
    filename : string
        Name of the file which should have the extension ``dst``.

    Returns
    -------
    data : tuple
        A tuple consisting of a list of taxa and a matrix.

    """

    try:
        f = open(filename)
    except:
        print("[!] Could not find the file {0}!".format(filename))

    taxa,matrix = [],[]
    
    
    for i,line in enumerate(f):
        if i > 0:
            taxa.append(line[0:10].strip())
            matrix.append([float(i) for i in
                re.split('\s+',line[11:].strip())])

    return taxa,matrix

