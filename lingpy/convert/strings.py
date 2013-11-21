# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-10-08 11:38
# modified : 2013-10-25 14:11
"""
Basic functions for the conversion of Python-internal data into strings.
"""

__author__="Johann-Mattis List"
__date__="2013-10-25"

import codecs
import regex as re
from ..settings import rcParams

def scorer2str(
        scorer
        ):
    """
    Convert a scoring function to a string.
    """
    
    # get sorted representation of characters
    chars = sorted(
            scorer.chars2int
            )

    # get the matrix
    matrix = scorer.matrix
    
    out = ''

    # write stuff to string
    for i,charA in enumerate(chars):
        out += charA
        for j,charB in enumerate(chars):
            out += '\t{0:.2f}'.format(scorer[charA,charB])
        out += '\n'

    return out

def msa2str(msa,wordlist=False):
    """
    Function converts an MSA object into a string.
    """
    
    out = ''

    # if wordlist ist set to True, don't write the header line and put the
    # after comment
    if wordlist: 
        formatter = "{0}\t{1:"+str(max([len(t) for t in msa['taxa']]))+'}'
        out += '# ' +msa['seq_id']+'\n'
        for a,b,c in zip(msa['ID'],msa['taxa'],msa['alignment']):
            out += formatter.format(a,b)+'\t'
            out += '\t'.join(c)+'\n'
        alm_len = len(c)

    elif type(msa) == dict:
        # get formatter
        formatter = '{0:'+str(max([len(t) for t in msa['taxa']]))+'}'
        out += msa['dataset']+'\n'
        out += msa['seq_id']+'\n'
        for a,b in zip(msa['taxa'],msa['alignment']):
            out += formatter.format(a)+'\t'
            out += '\t'.join(b)+'\n'
        alm_len = len(b)
    else:
        # get formatter
        formatter = '{0:'+str(max([len(t) for t in msa.taxa]))+'}'
        out += msa.dataset+'\n'
        out += msa.seq_id+'\n'
        for a,b in zip(msa.taxa,msa.alm_matrix):
            out += formatter.format(a)+'\t'
            out += '\t'.join(b)+'\n'
        alm_len = len(b)
    
    if 'local' in msa:
        local = msa['local']
    elif hasattr(msa,'local'):
        local = msa.local
    else:
        local = False

    if 'swaps' in msa:
        swaps = msa['swaps']
    elif hasattr(msa,'swaps'):
        swaps = msa.swaps
    else:
        swaps = False
    
    if local:
        if wordlist:
            out += formatter.format(0,"LOCAL")+'\t'
        else:
            out += formatter.format("LOCAL")+'\t'
        tmp = []
        for i in range(alm_len):
            if i in local:
                tmp += ['*']
            else:
                tmp += ['.']
        out += '\t'.join(tmp)+'\n'
    if swaps:
        if wordlist:
            out += formatter.format(0,"SWAPS")+'\t'
        else:
            out += formatter.format("SWAPS")+'\t'
        tmp = alm_len * ['.']
        for swap in swaps:
            a,b,c = swap
            tmp[a] = '+'
            tmp[b] = '-'
            tmp[c] = '+'
        out += '\t'.join(tmp)+'\n'

    return out

def matrix2dst(
        matrix,
        taxa = None,
        stamp = '',
        filename = '',
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
        if rcParams['verbose']: print(rcParams['fw'].format('dst')) # (filename,'dst'))

def pap2nex(
        taxa,
        paps,
        missing=0,
        filename=''
        ):
    """
    Function converts a list of paps into nexus file format.

    """
    out = '#NEXUS\n\nBEGIN DATA;\nDIMENSIONS ntax={0} NCHAR={1};\n'
    out += "FORMAT DATATYPE=STANDARD GAP=- MISSING={2} interleave=yes;\n"
    out += "MATRIX\n\n{3}\n;\n\nEND;"
    
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
    
    if not filename:
        return out.format(
                len(taxa),
                len(paps),
                missing,
                matrix
                )
    else:
        f = codecs.open(filename+'.nex','w')
        f.write(
                out.format(
                    len(taxa),
                    len(paps),
                    missing,
                    matrix
                    )
                )
        f.close()
        
        if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.nex'))
        return

def pap2csv(
        taxa,
        paps,
        filename=''
        ):
    """
    Write paps created by the Wordlist class to a csv-file.
    """

    out = "ID\t"+'\t'.join(taxa)+'\n'
    for key in sorted(paps,key=lambda x: int(re.sub(r'[^0-9]+','',str(x)))):
        out += '{0}\t{1}\n'.format(
            key,
            '\t'.join(str(i) for i in paps[key])
            )
    
    if not filename:
        return out
    else:
        f = codecs.open(filename+'.csv','w',"utf-8")
        f.write(out)
        f.close()

        if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.csv'))
        
        return
