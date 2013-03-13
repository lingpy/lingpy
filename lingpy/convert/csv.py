# created: Mo 21 Jan 2013 01:48:23  CET
# modified: Mo 21 Jan 2013 01:48:23  CET

__author__ = "Johann-Mattis List"
__date__ = "2013-01-21"

"""
Module provides functions and methods for the creation of csv-files.
"""

# imports
import re

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

    out = '# Wordlist\n'
    out += 'ID\t'+'\t'.join(header)+'\n'
    
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
    f.close()
    
    FileWriteMessage(filename,'csv').message('written')

    return 
