# created: Mo 21 Jan 2013 01:48:23  CET
# modified: Mo 21 Jan 2013 01:48:23  CET

__author__ = "Johann-Mattis List"
__date__ = "2013-01-21"

"""
Module provides functions and methods for the creation of csv-files.
"""

# imports
import re

def pap2csv(taxa,paps,filename='csv'):
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
    
    print("[i] Data has been written to file <{0}.csv>.".format(filename))
    
    return 
