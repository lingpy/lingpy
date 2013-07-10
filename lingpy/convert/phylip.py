# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-02 06:55
# modified : 2013-07-10 12:17
"""
Basic routines or creating Phylip output (distance matrices).
"""

__author__="Johann-Mattis List"
__date__="2013-07-10"


# external
import codecs

try:
    from ..algorithm.cython import misc
except:
    from ..algorithm.cython import _misc as misc

def matrix2dst(
        matrix,
        taxa = [],
        stamp = '',
        filename = '',
        verbose = True
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
    if not filename:
        return out
    else:
        f = codecs.open(filename+'.dst','w','utf-8')
        f.write(out)
        f.close()
        if verbose: print(FileWriteMessage(filename,'dst'))

