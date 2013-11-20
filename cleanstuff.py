# author   : Johann-Mattis List, Peter Bouda
# email    : mattis.list@gmail.com
# created  : 2013-07-10 14:22
# modified : 2013-07-20 12:26
"""
Script converts Python3-version of LingPy into a (hopefully valid) Python2-version.
"""

__author__="Johann-Mattis List, Peter Bouda"
__date__="2013-07-20"

#from glob import glob
import os
from os import path as osp
import codecs
import re


def cleanLP():
    parsD = []
    files = []
    
    for root, dirnames, filenames in os.walk('lingpy'):
        for f in filenames:
            if f.endswith('.py') or f.endswith('.pyx'):
                if 'settings.py' not in f:
                    files.append(os.path.join(root, f))
            elif len([x for x in ['.bin','.o','.so','.pyc','~','.swp'] if f.endswith(x)]) == 0:
                files.append(os.path.join(root, f))
    
    # create source target list for simple replacements
    st_list = [
            ('warning_empty_cons','W_empty_cons'),
            ('warning_failed_cons','W_failed_cons'),
            ('deprecation_warning','W_deprecation'),
            ('warning_deprecation','W_deprecation'),
            ('identical_scorer_warning','W_identical_scorer'),
            ('overwrite_scoring_function','W_overwrite_scorer'),
            ('warning_zero_division','E_zero_division'),
            ('warning_missing_module','W_missing_module'),
            ('empty_consensus_warning','W_empty_cons'),
            ('missing_module','W_missing_module'),
            ('sonority_consensus_warning','E_failed_cons'),
            ("rcParams['scale']","rcParams['align_scale']"),
            ("rcParams['factor']","rcParams['align_factor']"),
            ("rcParams['tree_calc']","rcParams['align_tree_calc']"),
            ("rcParams['gop']","rcParams['align_gop']"),
            ('rcParams["scale"]','rcParams["align_scale"]'),
            ('rcParams["factor"]','rcParams["align_factor"]'),
            ('rcParams["tree_calc"]','rcParams["align_tree_calc"]'),
            ('rcParams["gop"]','rcParams["align_gop"]'),
            ("rcParams['fw']","rcParams['M_file_written']"),
            ("rcParams['file_written']","rcParams['M_file_written']"),
            ]
    
    # iterate over each file and write a new version to the output
    for f in files:
        if not f.endswith('.bin'):
            stuff = codecs.open(f,'r','utf-8').read()
            oldstuff = stuff
            for source,target in st_list:
                stuff = stuff.replace(source,target)
            
            if stuff != oldstuff:
                print("[i] Converting file {0}...".format(f))
                nf = f.replace('lingpy','lingpy3')
                d = os.path.dirname(nf)
                if d and not os.path.exists(d):
                    os.makedirs(d)
                
                out = codecs.open(nf,'w','utf-8')
                out.write(stuff)
                out.close()
            
            oldstuff = stuff
            # find all rcParams and add them
            pars = re.findall(r'rcParams\[(.*?)\]',oldstuff)
            for p in pars:
                parsD += [p[1:-1]]
        else:
            pass
    
    return sorted(set(parsD))

if __name__ == "__main__":
    os.system('rm -r lingpy3/*')
    os.system('mkdir lingpy3')
    pars = cleanLP()

