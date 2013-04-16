# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-02 06:55
# modified : 2013-04-02 06:55
"""
Module provides functions and methods for the creation of csv-files.
"""

__author__="Johann-Mattis List"
__date__="2013-04-02"

# imports
import re
import json

from ..check.messages import FileWriteMessage
from .phylip import matrix2dst
from .misc import msa2str

def pap2csv(
        taxa,
        paps,
        filename='csv'
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

    f = open(filename+'.csv','w')
    f.write(out)
    f.close()

    FileWriteMessage(filename,'csv').message('written')
    
    return

def wl2csv(
        header,
        data,
        filename = 'csv',
        formatter = 'concept',
        **keywords
        ):
    """
    Write the basic data of a wordlist to file.
    """
    formatter = formatter.upper()

    # create output string
    out = '# Wordlist\n'

    # write meta to file
    if "meta" in keywords:
        meta = keywords["meta"]
    else:
        meta = {}
    
    kvpairs = {}
    jsonpairs = {}
    msapairs = {}
    distances = ''
    taxa = ''

    for k,v in meta.items():
        # simple key-value-pairs
        if type(v) in [str,int] or k == "tree":
            kvpairs[k] = v
        elif k == 'msa':
            for a,b in v.items():
                msapairs[a] = b
        elif k == 'distances':
            distances = matrix2dst(v,meta['taxa'])
        elif k == 'taxa':
            taxa = '\n'.join(meta['taxa'])
        else:
            jsonpairs[k] = v
    if kvpairs:
        out += '\n# META\n'
        for k,v in sorted(kvpairs.items(),key=lambda x:x[0]):
            out += '@{0}:{1}\n'.format(k,v)
    if taxa and keywords['taxa']:
        out += '\n# TAXA\n<taxa>\n'+taxa+'\n</taxa>\n'
    if jsonpairs:
        out += '\n# JSON\n'
        out += "<json>\n"
        out += json.dumps(jsonpairs,indent=4)
        out += '\n</json>\n'
    if msapairs:
        out += "\n# MSA\n"
        for k,v in msapairs.items():
            out += '#\n<msa id="{0}">\n'.format(k)
            out += msa2str(v)
            #out += v['seq_id']+'\n'
            #for t,alm in zip(v['taxa'],v['alignment']):
            #    out += t + '\t' + '\t'.join(alm)+'\n'
            out += "</msa>\n"
    if distances:
        out += '\n# DISTANCES\n<dst>\n'
        out += distances+'</dst>\n'

    out += '\n# DATA\nID\t'+'\t'.join(header)+'\n'
    
    # check for gloss in header to create nice output format
    if formatter in header:
        idx = header.index(formatter)
        formatter = None
    else:
        idx = False
        formatter = ''

    for key in sorted(data.keys()):
        
        # get the line
        line = data[key]
        
        # check for formatter
        if idx in range(len(line)):
            if line[idx] != formatter:
                out += '#\n'
                formatter = line[idx]

        # add the key 
        out += str(key)
        
        # add the rest of the values
        for value in line:
            if type(value) == list:
                try:
                    out += '\t'+' '.join(value)
                except:
                    out += '\t'+' '.join([str(v) for v in value])
            elif type(value) == int:
                out += '\t'+str(value)
            elif type(value) == float:
                out += '\t{0:.4f}'.format(value)
            else:
                out += '\t'+value
        out += '\n'

    f = open(filename + '.csv','w')
    f.write(out)
    if "stamp" in keywords:
        f.write(keywords['stamp'])
    f.close()
    FileWriteMessage(filename,'csv').message('written')

    return 
