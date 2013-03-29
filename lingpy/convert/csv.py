# created: Mo 21 Jan 2013 01:48:23  CET
# modified: Mo 21 Jan 2013 01:48:23  CET

__author__ = "Johann-Mattis List"
__date__ = "2013-01-21"

"""
Module provides functions and methods for the creation of csv-files.
"""

# imports
import re
import json

from ..check.messages import FileWriteMessage

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

    for k,v in meta.items():
        # simple key-value-pairs
        if type(v) in [str,int] or k == "tree":
            kvpairs[k] = v
            #out += '@{0}:{1}\n'.format(k,v)
        else:
            jsonpairs[k] = v
            #out += '\n# JSON\n<json id="{0}">\n'.format(k)
            #out += json.dumps(meta[k],indent=4)
            #out += "\n</json>\n\n"
    if kvpairs:
        out += '\n# META\n'
        for k,v in sorted(kvpairs.items(),key=lambda x:x[0]):
            out += '@{0}:{1}\n'.format(k,v)
    if jsonpairs:
        out += '\n# JSON\n'
        out += "<json>\n"
        out += json.dumps(jsonpairs,indent=4)
        out += '\n</json>\n'

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
