"""
Class for the handling of basic vocabulary lists.
"""
from __future__ import unicode_literals, print_function, division
from functools import partial
from glob import glob
import os

import numpy as np

from lingpy import log
from lingpy import util
from lingpy.basic.parser import QLCParserWithRowsAndCols
from lingpy.read.csv import csv2list


class BasVoc(QLCParserWithRowsAndCols):
    """
    Load a comparative collection of Swadesh lists (concepticon).

    Notes
    -----
    This collection may be useful for retrieving a subset of a given dataset,
    or for converting between conceptual items.

    Examples
    --------
    Load a BasVoc object without arguments in order to get the default object::

        >>> from lingpy.meaning import BasVoc
        >>> concepticon = BasVoc()

    Alternatively, load a pre-compiled object from LingPy::

        >>> from lingpy.meaning import concepticon

    Retrieve all original words in Jachontov's list concept list::

        >>> concepticon.get_list('jachontov','number','item')
        [['94', 'water'],
         ['25', 'eye'],
         ['45', 'know'],
         ['86', 'this'],
         ['84', 'tail'],
         ['87', 'thou'],
         ['28', 'fire'],
         ['89', 'tooth'],
         ['63', 'one'],
         ['32', 'full'],
         ['59', 'new'],
         ['42', 'I'],
         ['96', 'what'],
         ['82', 'sun'],
         ['61', 'nose'],
         ['37', 'hand'],
         ['18', 'dog'],
         ['24', 'egg'],
         ['81', 'stone'],
         ['88', 'tongue'],
         ['54', 'moon'],
         ['108', 'wind'],
         ['98', 'who'],
         ['104', 'salt'],
         ['50', 'louse'],
         ['91', 'two'],
         ['29', 'fish'],
         ['21', 'ear'],
         ['41', 'horn'],
         ['9', 'blood'],
         ['17', 'die'],
         ['110', 'year'],
         ['57', 'name'],
         ['10', 'bone'],
         ['33', 'give']]

    """
    def __init__(self, infile=None, col='list', row='key', conf=None):
        QLCParserWithRowsAndCols.__init__(
            self,
            infile or util.data_path('swadesh', 'swadesh.qlc'),
            row,
            col,
            conf or util.data_path('conf', 'swadesh.rc'))

        # get row and key index
        if not hasattr(self, '_rowidx'):
            # add indices to alias dictionary for swadesh lists
            for i, col in enumerate(self.cols):
                self._meta[col] = self._array[np.nonzero(self._array[:, i]), i][0]

    def get_dict(self, col='', row='', entry='', **keywords):
        """
        Return a dictionary representation for a given concept list.

        Parameters
        ----------
        col : str
            The concept list.
        row : str
            The concept (referenced by its unique ID).
        entry : str
            The entry that shall serve as key for the dictionary.

        Returns
        -------
        d : dict
            A dictionary with the unique IDs as key and the specified entry as
            value.
        """
        if row and not col:

            if row not in self.rows:
                raise ValueError("The row {0} you selected is not available!".format(row))
            else:
                data = self._dict[row]
                if not entry:
                    entries = data
                else:
                    entries = {}
                    for etr in entry.split(','):
                        idx = self._header[etr]

                        for key, value in data.items():
                            try:
                                entries[key] += [self[i][idx] for i in value]
                            except:
                                entries[key] = [self[i][idx] for i in value]

                return entries

        if col and not row:

            if col not in self.cols:
                raise ValueError("[!] The column you selected is not available!")
            else:
                data = {}
                for i, j in [
                    (self[i][self._rowIdx], i)
                    for i in self._array[:, self.cols.index(col)] if i != 0
                ]:
                    try:
                        data[i] += [j]
                    except:
                        data[i] = [j]

                if not entry:
                    entries = data
                else:
                    entries = {}
                    for etr in entry.split(','):
                        idx = self._header[etr]

                        for key, value in data.items():
                            try:
                                entries[key] += [self[i][idx] for i in value]
                            except KeyError:
                                entries[key] = [self[i][idx] for i in value]

                return entries

        elif row and col:
            raise ValueError("[!] You should specify only a value for row or for col!")
        else:
            for key in [k for k in keywords if k in self._alias]:
                if self._alias[key] == self._col_name:
                    entries = self.get_dict(col=keywords[key], entry=entry)
                    return entries

                elif self._alias[key] == self._row_name:
                    entries = self.get_dict(row=keywords[key], entry=entry)
                    return entries

            raise ValueError("[!] Neither rows nor columns are selected!")

    def get_list(self, swadlist, *entries):
        """
        Return a given concept list with the specified entries.

        Parameters
        ----------
        swadlist : str
            The concept list that shall be selected.
        *entries : str
            The entries that shall be included in the list.

        Returns
        -------
        l : list
            A list that contains the entries as specified.

        Examples
        --------

        >>> from lingpy.meaning import concepticon
        >>> lst = concepticon.get_list('jachontov', 'item')[:5]
        >>> lst
        ['water', 'eye', 'know', 'this', 'tail']

        """
        if entries:
            if len(entries) > 1:
                return [self[i, entries] for i in self._meta[swadlist]]
            else:
                return [self[i, entries[0]] for i in self._meta[swadlist]]
        else:
            return [self[i] for i in self._meta[swadlist]]

    def get_sublist(self, sublist, baselist, *entries):
        """
        Return the entries of one list that also occur in another list.

        Parameters
        ----------
        sublist : str
            The sublist whose entries shall be selected.
        baselst : str
            The name of the basic list of which the sublist shall be taken.
        *entries : str
            The entries ("item", "number", etc.) which shall be selected from
            the lists.

        Examples
        --------
        >>> from lingpy.meaning import concepticon
        >>> concepticon.get_sublist('dolgopolsky','jachontov','item')
        ['water',
         'eye',
         'thou',
         'tooth',
         'I',
         'what',
         'tongue',
         'who',
         'louse',
         'two',
         'name']

        Returns
        -------
        l : list
            A list containing the entries as specified.

        """
        listA = self.get_list(baselist)
        listB = self.get_list(sublist, 'key')

        if not entries:
            return [l for l in listA if l[self.header['key']] in listB]
        elif len(entries) > 1:
            return [
                [
                    l[self.header[x]] for x in entries
                ] for l in listA if l[self.header['key']] in listB
            ]
        else:
            return [l[self.header[entries[0]]]
                    for l in listA if l[self.header['key']] in listB]


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
                log.error(
                    "Multiple matches could be found: " + ', '.join(matches) + '...')
                raise
        else:
            return

    def __str__(self):
        return '\n'.join(self._keys)

    def __repr__(self):
        return "<concepticon with " + str(len(self._keys) - 1) + " conceptlists>"

    def compare(self, *lists, **keywords):
        """
        Compare multiple concept lists with each other.
        """
        kw = dict(mode='intersection', output=None, filename='output')
        kw.update(keywords)

        out = {}
        listlen = 0

        for k in lists:
            curlist = self[k]
            if not curlist:
                raise ValueError("Not all of the lists you selected are defined.")
            else:
                listlen += 1

        for i, k in enumerate(lists):
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
