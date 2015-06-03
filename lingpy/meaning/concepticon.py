# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-06-24 14:17
# modified : 2014-06-24 14:17
"""
This module provides an interface for a concepticon viewer.
"""

__author__="Johann-Mattis List"
__date__="2014-06-24"

from functools import partial

from ..settings import rcParams
from ..read.csv import csv2list
from glob import glob
import os
from .. import util
from .. import log


class Concepticon(object):

    def __init__(self):
        wlpath = partial(util.data_path, 'conceptlists')

        D = {}

        concepticon = csv2list(wlpath('concepticon.tsv'))
        header = [c.lower() for c in concepticon[0]][1:]
        D['base'] = {}
        D['base'][0] = header
        for line in concepticon[1:]:
            D['base'][line[0]] = line[1:]

        # create an artificial id for those entries which don't have one yet
        artid = -1

        for l in [lst for lst in glob(wlpath('*.tsv')) if 'concepticon' not in lst]:
            curlist = csv2list(l)

            # get the header
            header = [c.lower() for c in curlist[0]][:-1]

            name = os.path.split(l)[1][:-4].lower()

            D[name] = {}
            for line in curlist[1:]:
                if line[-1] == 'NA':
                    D[name][str(artid)] = line[:-1]
                    artid -= 1
                else:
                    D[name][line[-1]] = line[:-1]
            
            D[name][0] = header
    
        self._dict = D
        self._keys = list(self._dict.keys())
    
    def __getitem__(self, key):
        
        key = key.lower()

        # first, check whether key is a given concept list
        if key in self._keys:
            return self._dict[key]
        elif key.lower() in ' '.join(self._keys):
            matches = []
            for k in self._keys:
                if key in k: 
                    matches += [k]
            if len(matches) == 1:
                return self._dict[matches[0]]
            else:
                log.error("Multiple matches could be found: "+', '.join(matches)+'...')
                raise
        else:
            return 

    def __str__(self):
        
        return '\n'.join(self._keys)

    def __repr__(self):

        return "<concepticon with "+str(len(self._keys)-1)+" conceptlists>"
    
    def compare(self, *lists, **keywords):
        """
        Compare multiple concept lists with each other.
        """
        
        kw = dict(
                mode = 'intersection',
                output = None,
                filename = 'output',
                )
        kw.update(keywords)

        out = {}
        
        listlen = 0

        for k in lists:
            
            curlist = self[k]
            if not curlist:
                raise ValueError("Not all of the lists you selected are defined.")
            else:
                listlen += 1
        
        for i,k in enumerate(lists):
            
            curlist = self[k]
            header = curlist[0]
            gidx = header.index('gloss')
            for key in [x for x in curlist if x != 0]:
                gloss = curlist[key][gidx]
                try:
                    out[key][i] = gloss
                except KeyError:
                    out[key] = ['-' for k in range(listlen)]
                    out[key][i] = gloss

        for k in list(out.keys()):
            if kw['mode'] == 'intersection':
                if len([1 for x in out[k] if x != '-']) != listlen:
                    del out[k]
            elif kw['mode'] == 'difference':
                if len([1 for x in out[k] if x != '-']) == listlen:
                    del out[k]
            else:
                pass
            
        return out
