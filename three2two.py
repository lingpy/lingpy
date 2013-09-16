# author   : Johann-Mattis List, Peter Bouda
# email    : mattis.list@gmail.com
# created  : 2013-07-10 14:22
# modified : 2013-09-09 20:25
"""
Script converts Python3-version of LingPy into a (hopefully valid) Python2-version.
"""

__author__="Johann-Mattis List, Peter Bouda"
__date__="2013-09-09"

#from glob import glob
import os
from os import path as osp
import codecs
import re

def run3to2():
    files = []
    #inits = glob('lingpy/*')
    
    for root, dirnames, filenames in os.walk('lingpy'):
        for f in filenames:
            if f.endswith('.py') or f.endswith('.pyx'):
                files.append(os.path.join(root, f))
            elif len([x for x in ['.bin','.o','.so','.pyc','~','.swp'] if f.endswith(x)]) == 0:
                files.append(os.path.join(root, f))

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
            (' str:',' unicode:') 
            ]

    # pyxlist
    pyx_list = [
            (' str ',' unicode ')
            ]
    
    # iterate over each file and write a new version to the output
    for f in files:
        if not f.endswith('.bin'):
            print("[i] Converting file {0}...".format(f))
            stuff = codecs.open(f,'r','utf-8').read()
            for source,target in st_list:
                stuff = stuff.replace(source,target)

            # new specific line for str -> unicode
            stuff = re.sub(r'([\t \[])str([\t \n\(])',r'\1unicode\2',stuff)

            if f.endswith('pyx'):
                for source,target in pyx_list:
                    stuff = stuff.replace(source,target)
                stuff = re.sub("('.*?')",r"u\1",stuff)
                stuff = re.sub('(".*?")',r"u\1",stuff)
    
            nf = f.replace('lingpy','lingpy_build/lingpy') 
            d = os.path.dirname(nf)
            if d and not os.path.exists(d):
                os.makedirs(d)
            
            out = codecs.open(nf,'w','utf-8')
            if f.endswith('.py') or f.endswith('.pyx'):
                out.write(prefix)
            out.write(stuff)
            out.close()
        else:
            nf = f.replace('lingpy','lingpy_build/lingpy')
            d = os.path.dirname(nf)
            if d and not os.path.exists(d):
                os.makedirs(d)

            stuff = open(f,'rb').read()
            out = open(nf,'wb')
            out.write(stuff)
            out.close()

def run3to3():
    files = []

    for root, dirnames, filenames in os.walk('lingpy'):
        for f in filenames:
            if f.endswith('.py') or f.endswith('.pyx'):
                files.append(os.path.join(root, f))
            elif len([x for x in ['.bin','.o','.so','.pyc','~','.swp'] if f.endswith(x)]) == 0:
                files.append(os.path.join(root, f))
    
    for f in files:
        if not f.endswith('.bin'):
            print("[i] Copying file {0}...".format(f))
            stuff = codecs.open(f,'r','utf-8').read()
            
            nf = f.replace('lingpy','lingpy_build/lingpy') 
            d = os.path.dirname(nf)
            if d and not os.path.exists(d):
                os.makedirs(d)
            
            out = codecs.open(nf,'w','utf-8')
            out.write(stuff)
            out.close()
        else:
            nf = f.replace('lingpy','lingpy_build/lingpy')
            d = os.path.dirname(nf)
            if d and not os.path.exists(d):
                os.makedirs(d)

            stuff = open(f,'rb').read()
            out = open(nf,'wb')
            out.write(stuff)
            out.close()
