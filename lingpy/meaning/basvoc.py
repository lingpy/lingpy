# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-10-11 22:59
# modified : 2013-10-17 13:05
"""
Class for the handling of basic vocabulary lists.
"""

__author__="Johann-Mattis List"
__date__="2013-10-17"

from ..settings import rcParams
import os
from ..basic.parser import QLCParser
from six import text_type as str
import numpy as np

from .. import util


class BasVoc(QLCParser):
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
        QLCParser.__init__(
            self,
            infile or util.data_path('swadesh', 'swadesh.qlc'),
            conf or util.data_path('conf', 'swadesh.rc'))

        # get row and key index
        if not hasattr(self,'_rowidx'):
            try:
                rowIdx = self.header[self._alias[row]]
                colIdx = self.header[self._alias[col]]
            except:
                raise ValueError("[!] Could not find row and col in configuration or input file!")


            basic_rows = sorted(
                    set(
                        [self._data[k][rowIdx] for k in self._data if k != 0 and type(k) == int]
                        ),
                    key = lambda x: ('%s' % x).lower()
                    )
            basic_cols = sorted(
                    set(
                        [self._data[k][colIdx] for k in self._data if k != 0 and type(k) == int]
                        ),
                    key = lambda x: x.lower()
                    )
            
            # define rows and cols as attributes of the word list
            self.rows = basic_rows
            self.cols = basic_cols

            # define height and width of the word list
            self.height = len(self.rows)
            self.width = len(self.cols)
            
            # row and column index point to the place where the data of the main
            # items is stored in the original dictionary
            self._rowIdx = rowIdx
            self._colIdx = colIdx
            self._row_name = self._alias[row]
            self._col_name = self._alias[col]

            # create a basic array which assigns ids for the entries in a starling
            # manner. 

            # first, find out, how many items (== synonyms) are there maximally for
            # each row
            tmp_dict = {}
            for key,value in [(k,v) for k,v in self._data.items() if k != 0 and str(k).isnumeric()]:
                try:
                    tmp_dict[value[rowIdx]][value[colIdx]] += [key]
                except KeyError:
                    try:
                        tmp_dict[value[rowIdx]][value[colIdx]] = [key]
                    except KeyError:
                        tmp_dict[value[rowIdx]] = {}
                        tmp_dict[value[rowIdx]][value[colIdx]] = [key]
        
            # assign the values as _dict-attribute to the dictionary
            self._dict = tmp_dict

            # create the array by counting the maximal number of occurrences, store
            # the row names separately in a dictionary
            tmp_list = []
            row_dict = {}
            
            count = 0
            for k,d in self._dict.items():

                row_dict[k] = []

                # get maximal amount of "synonyms"
                m = max(
                        [len(x) for x in d.values()]
                        )

                for i in range(m):
                    tmp = []
                    for j in range(self.width):
                        try:
                            tmp.append(
                                    d[self.cols[j]][i]
                                    )
                        except:
                            tmp.append(0)
                    row_dict[k] += [count]
                    count += 1
                    tmp_list += [tmp]

            # create the array 
            self._array = np.array(tmp_list)
            self._idx = row_dict

            # add indices to alias dictionary for swadesh lists
            for i,col in enumerate(self.cols):
                self._meta[col] = self._array[np.nonzero(self._array[:,i]),i][0]
                
        # define a cache dictionary for stored data for quick access
        self._cache = {}
        
    def __getitem__(self,idx):
        """
        Method allows quick access to the data by passing the integer key.
        """
        try:
            return self._cache[idx]
        except:
            try:
                # return full data entry as list
                out = self._data[idx]
                self._cache[idx] = out
                return out
            except:
                try:
                    if type(idx[1]) in [tuple,list]:
                        out = [self._data[idx[0]][self._header[self._alias[i]]] for i in idx[1]]
                        self._cache[idx] = out
                        return out
                    else:
                        # return data entry with specified key word
                        out = self._data[idx[0]][self._header[self._alias[idx[1]]]]
                        return out
                except:
                    try:
                        if len(idx) == 2:
                            return [self._data[i][self._header[idx[1]]] for i in self._meta[idx[0]]]
                        else:
                            return [[self._data[i][self._header[x]] for x in idx[1:]] for i in self._meta[idx[0]]]
                    except:
                        raise KeyError("[ERROR] Key {0} cannot be interpreted.")

    def __getattr__(
            self,
            attr
            ):
        """
        Define how attributes are overloaded.
        """
        try:
            # get the right name
            nattr = self._alias[attr]

            if nattr == self._row_name:
                return self.rows
            elif nattr == self._col_name:
                return self.cols
            elif nattr in self._header:
                return self.get_entries(nattr)
            else:
                raise AttributeError("%r object has no attribute %r" %
                        (self.__class__,attr))
        except:
            try:
                return self._meta[attr] #return #return self._lists[attr]
            except:
                raise AttributeError("%r object has no attribute %r" %
                        (self.__class__,attr))

    def _clean_cache(self):

        """
        Function cleans the cache.
        """
        del self._cache
        self._cache = {}

    
    def get_dict(
            self,
            col = '',
            row = '',
            entry = '',
            **keywords
            ):
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
            try:
                return self._cache[row,entry]
            except:
                pass

            if row not in self.rows:
                print("[!] The row you selected is not available!")
            else:
                data = self._dict[row]
                if not entry:
                    entries = data
                else:
                    entries = {}
                    for etr in entry.split(','):
                        idx = self._header[etr]

                        for key,value in data.items():
                            try:
                                entries[key] += [self[i][idx] for i in value]
                            except:
                                entries[key] = [self[i][idx] for i in value]
                
                self._cache[row,entry] = entries
                return entries

        if col and not row:
            try:
                return self._cache[col,entry]
            except:
                pass

            if col not in self.cols:
                print("[!] The column you selected is not available!")
            else:
                data = {}
                for i,j in  [(self[i][self._rowIdx],i) for 
                        i in self._array[:,self.cols.index(col)] if i != 0]:
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

                        for key,value in data.items():
                            try:
                                entries[key] += [self[i][idx] for i in value]
                            except KeyError:
                                entries[key] = [self[i][idx] for i in value]
                
                self._cache[col,entry] = entries
                
                return entries
    
        elif row and col:
            print("[!] You should specify only a value for row or for col!")
        else:
            for key in [k for k in keywords if k in self._alias]:
                if self._alias[key] == self._col_name:
                    entries = self.get_dict(
                            col = keywords[key],
                            entry = entry,
                            )
                    self._cache[col,entry] = entries
                    return entries

                elif self._alias[key] == self._row_name:
                    entries = self.get_dict(
                            row = keywords[key],
                            entry = entry
                            )
                    self._cache[col,entry] = entries
                    return entries


            print("[!] Neither rows nor columns are selected!")
    def get_list(
            self,
            swadlist,
            *entries
            ):
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
                return [self[i,entries] for i in self._meta[swadlist]]
            else:
                return [self[i,entries[0]] for i in self._meta[swadlist]]
        else:
            return [self[i] for i in self._meta[swadlist]]
    
    def get_sublist(self,sublist,baselist,*entries):
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

        Example
        -------
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
        listB = self.get_list(sublist,'key')
        
        if not entries:
            return [l for l in listA if l[self.header['key']] in listB]
        elif len(entries) > 1:
            return [
                    [
                        l[self.header[x]] for x in entries
                        ] for l in listA if l[self.header['key']] in listB
                    ]
        else:
            return [l[self.header[entries[0]]] for l in listA if l[self.header['key']] in listB]

