# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-07-10 14:22
# modified : 2013-07-10 14:22
"""
Script converts Python3-version of LingPy into a (hopefully valid) Python2-version.
"""

__author__="Johann-Mattis List"
__date__="2013-07-10"

from glob import glob
import os
from os import path as osp
import codecs
import re

def run3to2():
    files = []
    inits = glob('lingpy/*')
    
    # read in all files
    while inits:
        f = inits.pop(0)
        if f.endswith('.py') or f.endswith('pyx'):
            files += [f]
        elif osp.isfile(f):
            #if f.endswith('.bin'):
            #    files += [f]
            if len([x for x in ['.bin','.o','.so','.pyc','~'] if f.endswith(x)]) == 0:
                files += [f]
        else:
            try:
                newfiles = glob(f+'/*')
                for f in newfiles:
                    inits += [f]
            except:
                pass
    
    prefix = """\
# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
"""
    
    # create source target list for simple replacements
    st_list = [
            ('input(','raw_input('),
            ]

    # pyxlist
    pyx_list = [
            (' str ',' unicode ')
            ]
    
    # iterate over each file and write a new version to the output
    for f in files:
        if not f.endswith('.bin'):
            stuff = codecs.open(f,'r','utf-8').read()
            for source,target in st_list:
                stuff = stuff.replace(source,target)
            if f.endswith('pyx'):
                for source,target in pyx_list:
                    stuff = stuff.replace(source,target)
                stuff = re.sub("('.*?')",r"u\1",stuff)
                stuff = re.sub('(".*?")',r"u\1",stuff)
    
            nf = f.replace('lingpy','lingpy2') #f.replace(path2package,'lingpy2')
    
            if '/' in nf:
                paths = nf.split('/')[:-1]
                for i in range(1,len(paths)+1):
                    if osp.isdir('/'.join(paths[:i])):
                        pass
                    else:
                        os.mkdir('/'.join(paths[:i]))
            
            out = codecs.open(nf,'w','utf-8')
            if f.endswith('.py') or f.endswith('.pyx'):
                out.write(prefix)
            out.write(stuff)
            out.close()
        else:
            nf = f.replace('lingpy','lingpy2')
            if '/' in nf:
                paths = nf.split('/')[:-1]
                for i in range(1,len(paths)+1):
                    if osp.isdir('/'.join(paths[:i])):
                        pass
                    else:
                        os.mkdir('/'.join(paths[:i]))

            stuff = open(f,'rb').read()
            out = open(nf,'wb')
            out.write(stuff)
            out.close()


