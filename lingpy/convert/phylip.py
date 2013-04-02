# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-02 06:55
# modified : 2013-04-02 06:55
"""
Basic routines or creating Phylip output (distance matrices).
"""

__author__="Johann-Mattis List"
__date__="2013-04-02"




try:
    from ..algorithm.cython import misc
except:
    from ..algorithm.cython import _misc as misc

def matrix2dst(
        matrix,
        taxa = False,
        stamp = False
        ):
    """
    Convert matrix to dst-format.
    """
    if not taxa:
        taxa = ['t_{0}'.format(i) for i in range(len(matrix))]

    out = ' {0}\n'.format(len(taxa))
    for i,taxon in enumerate(taxa):
        out += '{0:10}'.format(taxon)[0:11]
        out += ' '.join(['{0:2f}'.format(d) for d in
            matrix[i]])
        out += '\n'
    if stamp:
        out += '# {0}'.format(stamp)
    return out

#def wl2dst(
#        wl, # wordlist object
#        taxa = "taxa",
#        concepts = "concepts",
#        cognates = "cogid",
#        ):
#    """
#    Function converts wordlist to distance matrix.
#    """
#    # check for taxon attribute
#    taxa = getattr(wl,taxa)
#
#    # check for concepts
#    concepts = getattr(wl,concepts)
#
#    distances = []
#
#    for i,taxA in enumerate(taxa):
#        for j,taxB in enumerate(taxa):
#            if i < j:
#                
#                # get the two dictionaries
#                dictA = wl.get_dict(col=taxA,entry=cognates)
#                dictB = wl.get_dict(col=taxB,entry=cognates)
#    
#                # count amount of shared concepts
#                shared = 0
#                missing = 0
#                for concept in concepts:
#                    try:
#                        if [k for k in dictA[concept] if k in dictB[concept]]:
#                            shared += 1
#                        else:
#                            pass
#                    except KeyError:
#                        missing += 1
#    
#                # append to distances
#                distances += [ 1 - shared / (wl.height-missing)]
#    distances = misc.squareform(distances)
#    return distances
