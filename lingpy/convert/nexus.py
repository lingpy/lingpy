"""
Basic module for creating nexus output of linguistic data.
"""
__author__="Johann-Mattis List"
__date__="2012-12-05"

def pap2nex(taxa,paps,filename='nexus'):
    """
    Function converts a list of paps into nexus file format.


    """

    out = '#NEXUS\n\nBEGIN DATA;\nDIMENSIONS ntax={0} NCHAR={1};\n'
    out += "FORMAT DATATYPE=STANDARD GAP=- MISSING=? interleave=yes;\n"
    out += "MATRIX\n\n{2}\n;\n\nEND;"
    
    # get longest taxon
    maxTax = max([len(taxon) for taxon in taxa])

    # check whether paps are dict or list
    try:
        paps.keys()
        new_paps = []
        for key in paps:
            new_paps.append(paps[key])
    except:
        new_paps = paps


    # create the matrix
    matrix = ""
    
    for i,taxon in enumerate(taxa):
        tmp = '{0:XXX} '
        matrix += tmp.replace('XXX',str(maxTax)).format(taxon)
        matrix += ''.join([str(line[i]) for line in new_paps])
        matrix += '\n'

    f = open(filename+'.nex','w')
    f.write(
            out.format(
                len(taxa),
                len(paps),
                matrix
                )
            )

    print('[i] wrote stuff to file')
    
    return

        


