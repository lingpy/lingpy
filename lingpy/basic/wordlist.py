# *-* coding: utf-8 *-*
# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-14 00:21
# modified : 2014-12-02 21:10
"""
This module provides a basic class for the handling of word lists.
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

__author__="Johann-Mattis List"
__date__="2014-12-02"

import os
import numpy as np
import logging

from six import text_type as str

# basic lingpy imports
from ..read.qlc import read_qlc
from ..convert.strings import matrix2dst, pap2nex, pap2csv
from ..settings import rcParams
from .parser import QLCParser
from .ops import wl2dst, wl2dict, renumber, clean_taxnames, calculate_data, \
        wl2qlc, triple2tsv, tsv2triple

from ..algorithm import clustering as cluster
from ..algorithm import misc
from .. import util
from .. import log


def _write_file(filename, content, ext=None):
    if ext:
        filename = filename + '.' + ext
    util.write_text_file(filename, content)


class Wordlist(QLCParser):
    """
    Basic class for the handling of multilingual word lists.

    Parameters
    ----------
    filename : { string, dict }
        The input file that contains the data. Otherwise a dictionary with
        consecutive integers as keys and lists as values with the key 0
        specifying the header.

    row : str (default = "concept")
        A string indicating the name of the row that shall be taken as the
        basis for the tabular representation of the word list.
    
    col : str (default = "doculect")
        A string indicating the name of the column that shall be taken as the
        basis for the tabular representation of the word list.
    
    conf : string (default='')
        A string defining the path to the configuration file. 

    Notes
    -----
    A word list is created from a dictionary containing the data. Two keywords
    (row and col) define, which of the dimensions of the original data should
    be used as row and as column of the tabular display. A configuration file
    can be used to change basic names and aliases for the data being used, and
    the classes (data types) of the entries.

    A couple of methods is provided along with the word list class in order to
    access the multi-dimensional input data. The main idea is to provide an
    easy way to access two-dimensional slices of the data by specifying which
    entry type should be returned. Thus, if a word list consists not only of
    simple orthographical entries but also of IPA encoded phonetic
    transcriptions, both the orthographical source and the IPA transcriptions
    can be easily accessed as two separate two-dimensional lists.
    
    """

    def __init__(self, filename, row='concept', col='doculect', conf=None):
        
        QLCParser.__init__(self, filename, conf or util.data_path('conf', 'wordlist.rc'))
        
        # check whether additional data has to be loaded
        if not hasattr(self,'_rowidx'):
            # check whether additional data has to be loaded
            # retrieve basic types for rows and columns from the word list
            try:
                rowIdx = self.header[self._alias[row]]
                colIdx = self.header[self._alias[col]]
            except:
                raise ValueError("[!] Could not find row and col in configuration or input file!")
            
            basic_rows = sorted(
                    set(
                        [self._data[k][rowIdx] for k in self._data if k != 0 and isinstance(k, int)]
                        ),
                    key = lambda x: x.lower()
                    )
            basic_cols = sorted(
                    set(
                        [self._data[k][colIdx] for k in self._data if k != 0 and isinstance(k, int)]
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
                
        # define a cache dictionary for stored data for quick access
        self._cache = {}

        # check for taxa in meta
        if 'taxa' in self._alias:
            if self._alias['taxa'] not in self._meta:
                self._meta[self._alias['taxa']] = self.cols

    def __getitem__(
            self,
            idx
            ):
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
                    # return data entry with specified key word
                    out = self._data[idx[0]][self._header[self._alias[idx[1]]]]
                    self._cache[idx] = out
                    return out
                except:
                    try:
                        out = self._meta[idx]
                        self._cache[idx] = out
                        return out
                    except:
                        pass

    def _clean_cache(self):

        """
        Function cleans the cache.
        """
        del self._cache
        self._cache = {}

    def add_entries(
            self,
            entry,
            source,
            function,
            override=False,
            **keywords
            ):
        """
        Add new entry-types to the word list by modifying given ones.

        Parameters
        ----------
        entry : string
            A string specifying the name of the new entry-type to be added to the
            word list.

        source : string
            A string specifying the basic entry-type that shall be modified. If
            multiple entry-types shall be used to create a new entry, they
            should be passed in a simple string separated by a comma.

        function : function
            A function which is used to convert the source into the target
            value.

        keywords : {dict}
            A dictionary of keywords that are passed as parameters to the
            function.

        Notes
        -----
        This method can be used to add new entry-types to the data by
        converting given ones. There are a lot of possibilities for adding new
        entries, but the most basic procedure is to use an existing entry-type
        and to modify it with help of a function.

        """
        self._add_entries(entry, source, function, override, **keywords)
        self._clean_cache()

    def get_dict(
            self,
            col='',
            row='',
            entry='',
            **keywords
            ):
        """
        Function returns dictionaries of the cells matched by the indices.

        Parameters
        ----------
        col : string (default="")
            The column index evaluated by the method. It should contain one of the
            values in the row of the
            :py:class:`~lingpy.basic.wordlist.Wordlist` instance, usually a
            taxon (language) name.
        
        row : string (default="")
            The row index evaluated by the method. It should contain one of the
            values in the row of the
            :py:class:`~lingpy.basic.wordlist.Wordlist` instance, usually a
            taxon (language) name.
        
        entry : string (default="")
            The index for the entry evaluated by the method. It can be used to specify
            the datatype of the rows or columns selected. As a default, the
            indices of the entries are returned.

        Returns
        -------
        entries : dict
            A dictionary of keys and values specifying the selected part of the
            data. Typically, this can be a dictionary of a given language with
            keys for the concept and values as specified in the "entry"
            keyword.

        Notes
        -----
        The "col" and "row" keywords in the function are all aliased according
        to the description in the ``wordlist.rc`` file. Thus, instead of
        using these attributes, the aliases can also be taken. For selecting a
        language, one may type something like::

            >>> Wordlist.get_dict(language='LANGUAGE')
        
        and for the selection of a concept, one may type something like::

            >>> Wordlist.get_dict(concept='CONCEPT')    
        
        See the examples below for details.

        Examples
        --------
        Load the ``harry_potter.csv`` file::
            
            >>> wl = Wordlist('harry_potter.csv')

        Select all IPA-entries for the language "German"::
        
            >>> wl.get_dict(language='German',entry='ipa')
            {'Harry': ['haralt'], 'hand': ['hant'], 'leg': ['bain']}
        
        Select all words (orthographical representation) for the concept
        "Harry"::
        
            >>> wl.get_dict(concept="Harry",entry="words")
            {'English': ['hæri'], 'German': ['haralt'], 'Russian': ['gari'], 'Ukrainian': ['gari']}
            
        Note that the values of the dictionary that is returned are always
        lists, since it is possible that the original file contains synonyms
        (multiple words corresponding to the same concept).

        See also
        --------
        Wordlist.get_list
        Wordlist.add_entries

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
                    idx = self._header[entry]

                    for key,value in data.items():
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
                    idx = self._header[entry]

                    for key,value in data.items():
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
            row='',
            col='',
            entry='',
            flat=False,
            **keywords
            ):
        """
        Function returns lists of rows and columns specified by their name.

        Parameters
        ----------
        row: string (default = '')
            The row name whose entries are selected from the data.

        col : string (default = '')
            The column name whose entries are selected from the data.
        
        entry: string (default = '')
            The entry-type which is selected from the data.

        flat : bool (default = False)
            Specify whether the returned list should be one- or
            two-dimensional, or whether it should contain gaps or not.

        Returns
        -------
        data : list
            A list representing the selected part of the
            data.


        Notes
        -----
        The 'col' and 'row' keywords in the function are all aliased according
        to the description in the ``wordlist.rc`` file. Thus, instead of
        using these attributes, the aliases can also be taken. For selecting a
        language, one may type something like::

            >>> Wordlist.get_list(language='LANGUAGE')
        
        and for the selection of a concept, one may type something like::

            >>> Wordlist.get_list(concept='CONCEPT')    
        
        See the examples below for details.
        
        Examples
        --------
        Load the ``harry_potter.csv`` file::
            
            >>> wl = Wordlist('harry_potter.csv')

        Select all IPA-entries for the language "German"::
        
            >>> wl.get_list(language='German',entry='ipa'
            ['bain', 'hant', 'haralt']

        Note that this function returns 0 for missing values (concepts that
        don't have a word in the given language). If one wants to avoid this,
        the 'flat' keyword should be set to *True*.
        
        Select all words (orthographical representation) for the concept
        "Harry"::
        
            >>> wl.get_list(concept="Harry",entry="words")
            [['Harry', 'Harald', 'Гари', 'Гарi']]            
        
        Note that the values of the list that is returned are always
        two-dimensional
        lists, since it is possible that the original file contains synonyms
        (multiple words corresponding to the same concept). If one wants to
        have a flat representation of the entries, the 'flat' keyword should be
        set to *True*::

            >>> wl.get_list(concept="Harry",entry="words",flat=True)
            ['hæri', 'haralt', 'gari', 'hari']
            
        See also
        --------
        Wordlist.get_list
        Wordlist.add_entries

        """
        # if row is chosen
        if row and not col:
            # first try to return what's in the cache
            try:
                return self._cache[row,entry,flat]
            except:
                pass

            # otherwise, start searching
            if row not in self.rows:
                print("[!] The row you selected is not available.")
            else:
                # first, get the row ids
                data = self._array[self._idx[row]]

                # if only row is chosen, return the ids
                if not entry:

                    # check for flat representation
                    if flat:
                        entries = [i for i in data.flatten() if i != 0]
                    else:
                        entries = data.tolist()

                # if row and entry-type is chosen, return the entry-type
                else:
                    # get the index for the entry in the data dictionary
                    idx = self._header[entry]
                    
                    if flat:
                        # get the entries
                        entries = [self[i][idx] for i in 
                                data.flatten() if i != 0]
                    else:
                        # get the entries
                        entries = data.tolist()
                        for i,line in enumerate(data):
                            for j,cell in enumerate(line):
                                if cell != 0:
                                    entries[i][j] = self[cell][idx]
                    
                # append entries to cache
                self._cache[row,entry,flat] = entries

                return entries

        # if column is chosen
        elif col and not row:
            # try to return the cache
            try:
                return self._cache[col,entry,flat]
            except:
                pass

            if col not in self.cols:
                print("[!] The column you selected is not available!")
                return
            else:
                data = self._array[:,self.cols.index(col)]
                
                if not entry:
                    if flat:
                        entries = [i for i in data if i != 0]
                    else:
                        entries = data.tolist()
                else:
                    idx = self._header[entry]

                    if flat:
                        entries = [self[i][idx] for i in data if i != 0]

                    else:
                        entries = []
                        for i in data:
                            if i != 0:
                                entries.append(self[i][idx])
                            else:
                                entries.append(0)
            
            self._cache[col,entry,flat] = entries
            
            return entries
        
        elif row and col:
            print("[!] You should specify only a value for row or for col!")
        else:
            for key in [k for k in keywords if k in self._alias]:
                if self._alias[key] == self._col_name:
                    entries = self.get_list(
                            col = keywords[key],
                            entry = entry,
                            flat = flat
                            )
                    self._cache[col,entry,flat] = entries
                    return entries

                elif self._alias[key] == self._row_name:
                    entries = self.get_list(
                            row = keywords[key],
                            entry = entry,
                            flat = flat
                            )
                    self._cache[col,entry,flat] = entries
                    return entries

            print("[!] Neither rows nor columns are selected!")
            return 
    
    def get_entries(
            self,
            entry
            ):
        """
        Return all entries matching the given entry-type as a two-dimensional list.

        Parameters
        ----------
        entry : string
            The entry-type of the data that shall be returned in tabular
            format.

        """

        try:
            return self._cache[self._alias[entry]]
        except:
            pass

        if entry in self._header:
            
            # get the index
            idx = self._header[entry]

            entries = []

            for row in self._array:
                tmp = [0 for i in row]
                for i,cell in enumerate(row):
                    if cell != 0:
                        tmp[i] = self[cell][idx]
                entries.append(tmp)

            # add entries to cache
            self._cache[self._alias[entry]] = entries

            return entries

    def get_etymdict(
            self,
            ref="cogid",
            entry='',
            loans=False,
            fuzzy=False
            ):
        """
        Return an etymological dictionary representation of the word list.

        Parameters
        ----------
        ref : string (default = "cogid")
            The reference entry which is used to store the cognate ids.

        entry : string (default = '')
            The entry-type which shall be selected.

        loans : bool (default=False)
            Treat all negative cognate ids as if they were positive ones.
            Negative IDs can be used to indicate borrowings, or, more
            accurately, xenologs.

        fuzzy : bool (default=False)
            d

        Returns
        -------
        etym_dict : dict
            An etymological dictionary representation of the data.

        Notes
        -----
        In contrast to the word-list representation of the data, an
        etymological dictionary representation sorts the counterparts according to
        the cognate sets of which they are reflexes. 

        """
        
        # make an alias for the reference
        ref = (self._alias[ref],loans)

        # make converting function for loans
        if loans:
            f = lambda x: abs(x)
        else:
            f = lambda x: x

        # check in the cache
        try:
            return self._cache[ref,entry]
        except:
            pass

        # create an etymdict object
        try:
            self._etym_dict[ref]
        except:
            try:
                self._etym_dict
                self._etym_dict[ref] = {}
            except:
                self._etym_dict = {ref:{}}
            
            # get the index for the cognate id 
            cogIdx = self._header[ref[0]]
            
            # iterate over all data
            if not fuzzy:
                for key in self:
                    cogid = f(self[key][cogIdx])
                    colIdx = self.cols.index(self[key][self._colIdx])

                    # assign new line if cogid is missing
                    if cogid not in self._etym_dict[ref]:
                        self._etym_dict[ref][cogid] = [0 for i in range(self.width)]
                    
                    # assign the values for the current session
                    try:
                        self._etym_dict[ref][cogid][colIdx] += [key]
                    except:
                        self._etym_dict[ref][cogid][colIdx] = [key]

            elif fuzzy:
                for key in self:
                    cogids = self[key][cogIdx]
                    colIdx = self.cols.index(self[key][self._colIdx])
                    for cog in cogids:
                        cogid = f(cog)
                        if cog not in self._etym_dict[ref]:
                            self._etym_dict[ref][cog] = [0 for i in range(self.width)]

                        # assign new values for the current session
                        try:
                            self._etym_dict[ref][cog][colIdx] += [key]
                        except:
                            self._etym_dict[ref][cog][colIdx] = [key]
                        
        
        if entry:
            # create the output
            etym_dict = {}

            # get the index of the header
            idx = self._header[entry]

            # retrieve the values 
            for key,values in self._etym_dict[ref].items():
                etym_dict[key] = []
                for value in values:
                    if value != 0:
                        etym_dict[key].append(
                                [self[v][idx] for v in value]
                                )
                    else:
                        etym_dict[key].append(0)
        else:
            etym_dict = self._etym_dict[ref]
        
        # add the stuff to the cache
        self._cache[ref,entry] = etym_dict

        return etym_dict

    def get_paps(
            self,
            ref='cogid',
            entry='concept',
            missing=0,
            loans=False
            ):
        """
        Function returns a list of present-absent-patterns of a given word list.

        Parameters
        ----------
        ref : string (default = "cogid")
            The reference entry which is used to store the cognate ids.
        entry : string (default = "concept")
            The field which is used to check for missing data.
        missing : string,int (default = 0)
            The marker for missing items.
        """
         
        try:
            return self._cache['#paps#'+str(missing)+'#',ref]
        except:
            pass
        
        etym_dict = self.get_etymdict(ref=ref, entry=entry, loans=loans)            

        # create dictionary for paps
        paps = {}

        # create dictionary that stores missing data
        missed = {}

        # retrieve the values
        for key,values in etym_dict.items():
            paps[key] = []

            # check for missing data
            meanings = set()
            for value in values:
                if value:
                    for v in value:
                        meanings.add(v)
            if len(meanings) == 1:
                meaning = meanings.pop()

                if meaning not in missed:
                    
                    # get the list in the wordlist of self
                    tmp = np.array(self.get_list(row=meaning))
                    
                    # get the sum of the list
                    tmp = sum(tmp)
                    
                    # get all languages which are zero
                    gaps = [i for i in range(self.width) if not tmp[i]]

                    # append gaps to missing
                    missed[meaning] = gaps
            else:
                meaning = False

            for i,value in enumerate(values):
                if value:
                    paps[key].append(1)
                else:
                    if meaning:
                        if i in missed[meaning]:
                            paps[key].append(missing)
                        else:
                            paps[key].append(0)
                    else:
                        paps[key].append(1)
        
        self._cache['#paps#'+str(missing)+'#',ref] = paps
        
        return paps
    
    def calculate(
            self,
            data,
            taxa='taxa',
            concepts='concepts',
            ref='cogid',
            **keywords
            ):
        """
        Function calculates specific data.

        Parameters
        ----------
        data : str
            The type of data that shall be calculated. Currently supports

            * "tree": calculate a reference tree based on shared cognates
            * "dst": get distances between taxa based on shared cognates
            * "cluster": cluster the taxa into groups using different methods

        """
        
        calculate_data(self,data,taxa,concepts,ref,**keywords)

    def renumber(
            self,
            source,
            target='',
            override=False
            ):
        """
        Renumber a given set of string identifiers by replacing the ids by integers.
        
        Parameters
        ----------
        source : str
            The source column to be manipulated.

        target : str (default='')
            The name of the target colummn. If no name is chosen, the target
            column will be manipulated by adding "id" to the name of the source
            column.

        override : bool (default=False)
            Force to overwrite the data if the target column already exists.

        Notes
        -----
        In addition to a new column, an further entry is added to the "_meta"
        attribute of the wordlist by which newly coined ids can be retrieved
        from the former string attributes. This attribute is called
        "source2target" and can be accessed either via the "_meta" dictionary
        or directly as an attribute of the wordlist.

        """
        renumber(self,source,target,override=override)

    def _clean(
            self,
            source,
            f=lambda x: ''.join([t for t in x if t not in '()[] {},;:'])
            ):
        """
        Clean given names of doculects to make them work in Newick.

        Notes
        -----
        The function basically removes all brackes, whitespace, and the like
        from the taxon names in a wordlist. This is quite important for later
        calculation of trees and the like, since the Newick representation
        format requires taxon-names to contain no brackets, and other
        characters such as colons or dots. It is also useful for display, since
        language names like "German (Standard)" do not look really appealing in
        basic applications.
        """
        
        if self._alias[source] == 'doculect':
            clean_taxnames(self, self._alias[source], f)
    
    def _output(
            self,
            fileformat,
            **keywords
            ):
        """
        Internal function that eases its modification by daughter classes.
        """

        # check for stamp attribute
        if hasattr(self,"_stamp"):
            keywords["stamp"] = self._stamp
        else:
            keywords["stamp"] = ''

        # add the default parameters, they will be checked against the keywords
        defaults = {
                'cols'       : False,
                'distances'  : False,
                'entries'    : ("concept","counterpart"),
                'entry'      : 'concept',
                'fileformat' : fileformat,
                'filename'   : rcParams['filename'],
                'formatter'  : 'concept',
                'loans'      : False,
                'meta'       : self._meta,
                'missing'    : 0,
                'ref'        : 'cogid',
                'rows'       : False,
                'subset'     : False, # setup a subset of the data,
                'taxa'       : 'taxa',
                'threshold'  : 0.6, # threshold for flat clustering
                'tree_calc'  : 'neighbor',
                }
            
        # compare with keywords and add missing ones
        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        if fileformat in ['triple','triples','triples.tsv']:
            tsv2triple(self, keywords['filename']+'.'+fileformat)

        elif fileformat in ['paps.nex','paps.csv']:
            paps = self.get_paps(
                    ref=keywords['ref'],
                    entry=keywords['entry'],
                    missing=keywords['missing']
                    )
            if fileformat == 'paps.nex':
                pap2nex(
                        self.taxa,
                        paps,
                        missing=keywords['missing'],
                        filename=keywords['filename']+'.paps'
                        )
            elif fileformat == 'paps.csv':
                pap2csv(
                        self.cols,
                        paps,
                        filename=keywords['filename']+'.paps'
                        )
        
        # simple printing of taxa
        elif fileformat == 'taxa':
            if hasattr(self,'taxa'):
                out = ''
                for taxon in self.taxa:
                    out += taxon + '\n'
                util.write_text_file(keywords['filename'] + '.taxa', out)
            else:
                raise ValueError(
                        "[i] Taxa are not available."
                        )
        
        # csv-output
        elif fileformat in ['csv','qlc','tsv']:
            if fileformat == 'csv':
                log.deprecated('csv','qlc')

            # get the header line
            header = sorted(
                    [s for s in set(self._alias.values()) if s in self._header],
                    key = lambda x: self._header[x]
                    )
            header = [h.upper() for h in header]

            # check for taxa in meta
            if not 'taxa' in self._meta:
                self._meta['taxa'] = self.taxa

            # get the data, in case a subset is chosen
            if not keywords['subset']:
                # write stuff to file
                wl2qlc(
                        header,
                        self._data,
                        **keywords
                        )
            else:
                cols,rows = keywords['cols'],keywords['rows']

                if not isinstance(cols, (list, tuple, bool)):
                    raise ValueError(
                            "[i] Argument 'cols' should be list or tuple."
                            )
                if not isinstance(rows, (dict,bool)):
                    raise ValueError(
                            "[i] Argument 'rows' should be a dictionary."
                            )

                # check for chosen header
                if cols:
                    # get indices for header
                    indices = [self._header[x] for x in cols]
                    header = [c.upper() for c in cols]
                else:
                    indices = [r for r in range(len(self.header))]

                if rows:
                    stmts = []
                    for key,value in rows.items():
                        if key == 'ID':
                            stmts += ["key "+value] #
                        else:
                            idx = self._header[key]
                            stmts += ["line[{0}] ".format(idx)+value]

                log.debug("calculated what should be excluded")

                # get the data
                out = {}

                for key,line in self._data.items():
                    if log.get_level() <= logging.DEBUG:
                        print(key)

                    if rows:
                        if eval(" and ".join(stmts)):
                            out[key] = [line[i] for i in indices]
                    else:
                        out[key] = [line[i] for i in indices]
                
                log.debug("passing data to wl2qlc")
                wl2qlc(
                        header,
                        out,
                        **keywords
                        )
        
        # output dst-format (phylip)
        elif fileformat == 'dst':
            # check for distances as keyword
            if 'distances' not in self._meta:
                self._meta['distances'] = wl2dst(self,**keywords)
            
            # write data to file
            out = matrix2dst(
                    self._meta['distances'],
                    self.taxa,
                    stamp=keywords['stamp']
                    )
            _write_file(keywords['filename'], out, fileformat)
        # output tre-format (newick)
        elif fileformat in ['tre','nwk']: #,'cluster','groups']:
            if 'tree' not in self._meta:
            
                # check for distances
                if 'distances' not in self._meta:
                    self._meta['distances'] = wl2dst(self)
                    
                if keywords['tree_calc'] == 'neighbor':
                    tree = cluster.neighbor(
                            self._meta['distances'],
                            self.cols,
                            distances=keywords['distances']
                            )
                elif keywords['tree_calc'] == 'upgma':
                    tree = cluster.ugpma(
                            self._meta['distances'],
                            self.cols,
                            distances=keywords['distances']
                            )
            else:
                tree = self._meta['tree']

            _write_file(keywords['filename'], '{0}'.format(tree), fileformat)

        elif fileformat in ['cluster','groups']:
            if 'distances' not in self._meta:
                
                self._meta['distances'] = wl2dst(self) # check for keywords

            if 'groups' not in self._meta:
                self._meta['groups'] = cluster.matrix2groups(
                        keywords['threshold'],
                        self._meta['distances'],
                        self.taxa
                        )
            lines = []
            for taxon, group in sorted(self._meta['groups'].items(), key=lambda x: x[0]):
                lines.append('{0}\t{1}'.format(taxon, group))
            _write_file(keywords['filename'], '\n'.join(lines), fileformat)

        elif fileformat in ['starling','star.csv']:
            # make lambda inline for data-check
            l = lambda x: ['-' if x == 0 else x][0]

            lines = []
            if 'cognates' not in keywords:
                lines.append('ID\tConcept\t'+'\t'.join(self.taxa))
                for i, concept in enumerate(self.concepts):
                    for line in self.get_list(row=concept,entry=keywords['entry']):
                        lines.append(
                            str(i+1)+'\t'+concept+'\t'+'\t'.join([l(t) for t in line]))
            else:
                lines.append(
                    'ID\tConcept\t'+'\t'.join(['{0}\t COG'.format(t) for t in self.taxa]))
                for i, concept in enumerate(self.concepts):
                    cogs = self.get_list(row=concept,entry=keywords['cognates'])
                    for j,line in enumerate(self.get_list(row=concept,entry=keywords['entry'])):
                        part = '\t'.join('{0}\t{1}'.format(l(a),b) for a,b in zip(line,cogs[j]))
                        lines.append(str(i+1)+'\t'+concept+'\t'+part)

            _write_file(
                keywords['filename'],
                '\n'.join(lines),
                'starling_' + keywords['entry'] + '.csv')

        elif fileformat == 'separated':
            if not os.path.isdir(keywords['filename']):
                os.mkdir(keywords['filename'])

            for l in self.cols:
                if 'ignore_keys' in keywords:
                    lines = ['']
                else:
                    lines = ['ID\t']
                lines[0] += '\t'.join([x.upper() for x in keywords['entries']])
                for key in self.get_list(col=l,flat=True):
                    line = [key] if 'ignore_keys' not in keywords else []
                    for entry in keywords['entries']:
                        tmp = self[key,entry]
                        if isinstance(tmp, list):
                            tmp = ' '.join([str(x) for x in tmp])
                        
                        line += [tmp]
                    lines.append('\t'.join('{0}'.format(x) for x in line))
                _write_file(
                    '{0}/{1}'.format(keywords['filename'], l),
                    '\n'.join(lines),
                    'tsv')

    def output(
            self,
            fileformat,
            **keywords
            ):
        """
        Write wordlist to file.

        Parameters
        ----------
        fileformat : {"qlc", "tre","nwk","dst", "taxa", "starling", "paps.nex", "paps.csv"}
            The format that is written to file. This corresponds to the file
            extension, thus 'qlc' creates a file in qlc-format, 'dst' creates
            a file in Phylip-distance format, etc.
        filename : str
            Specify the name of the output file (defaults to a filename that
            indicates the creation date).
        subset : bool (default=False)
            If set to c{True}, return only a subset of the data. Which subset
            is specified in the keywords 'cols' and 'rows'.
        cols : list
            If *subset* is set to c{True}, specify the columns that shall be
            written to the csv-file.
        rows : dict
            If *subset* is set to c{True}, use a dictionary consisting of keys
            that specify a column and values that give a Python-statement in
            raw text, such as, e.g., "== 'hand'". The content of the specified
            column will then be checked against statement passed in the
            dictionary, and if it is evaluated to c{True}, the respective row
            will be written to file.
        cognates : str
            Name of the column that contains the cognate IDs if 'starling' is
            chosen as an output format.

        missing : { str, int } (default=0)
            If 'paps.nex' or 'paps.csv' is chosen as fileformat, this character
            will be inserted as an indicator of missing data.

        tree_calc : {'neighbor', 'upgma'}
            If no tree has been calculated and 'tre' or 'nwk' is chosen as
            output format, the method that is used to calculate the tree.

        threshold : float (default=0.6)
            The threshold that is used to carry out a flat cluster analysis if
            'groups' or 'cluster' is chosen as output format.

        """
        return self._output(fileformat,**keywords)
    
    def _export(
            self,
            fileformat,
            sections=None,
            entries=None,
            entry_sep='',
            item_sep='',
            template='',
            exclude=None,
            entry_start='',
            entry_close='',
            **keywords
            ):
        """
        Export a wordlist to various file formats.
        """
        # check for sections
        if not sections:
            if fileformat == 'txt':
                sections = dict(
                        h1 = ('concept', '\n# Concept: {0}\n'),
                        h2 = ('cogid', '## Cognate-ID: {0}\n'),
                        )
            elif fileformat == 'tex':
                sections = dict(
                        h1 = ('concept', r'\section{{Concept: ``{0}"}}'+'\n'),
                        h2 = ('cogid', r'\subsection{{Cognate Set: ``{0}"}}'+'\n')
                        )
            elif fileformat == 'html':
                
                sections = dict(
                        h1 = ('concept', '<h1>Concept: {0}</h1>'),
                        h2 = ('cogid', '<h2>Cognate Set: {0}</h2>')
                        )
    
        # check for entries
        if not entries:
            if fileformat == 'txt':
                entries = [
                        ('language', '{0} '),
                        ('ipa', '{0}\n')
                        ]
            elif fileformat == 'tex':
                entries = [
                        ('language', '{0} '),
                        ('ipa', '[{0}]'+'\n')
                        ]
            elif fileformat == 'html':
                entries = [
                        ('language', '{0}&nbsp;'),
                        ('ipa', '[{0}]\n')
                        ]
        
        # setup defaults
        defaults = dict(
                filename = rcParams['filename'] 
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        # get the temporary dictionary
        out = wl2dict(
                self,
                sections,
                entries,
                exclude
                )
    
        # assign the output string
        out_string = ''
    
        # iterate over the dictionary and start to fill the string
        for key in sorted(out, key=lambda x: str(x).lower()):
            
            # write key to file
            out_string += key[1]
            
            # reassign tmp
            tmp = out[key]
    
            # set the pointer and the index
            pointer = {0:[tmp,sorted(tmp.keys())]}
            

            break_loop = False

            while True:
                
                if break_loop:
                    break

                idx = max(pointer.keys())

                # check for type of current point
                if isinstance(tmp, dict):
                    
                    if pointer[idx][1]:
                        next_key = pointer[idx][1].pop()
                        out_string += next_key[1]
                        tmp = pointer[idx][0][next_key]
                        if isinstance(tmp, dict):
                            pointer[idx+1] = [tmp,sorted(tmp.keys())]
                        else:
                            pointer[idx+1] = [tmp,tmp]
                    else:
                        del pointer[idx]
                        if idx == 0:
                            break_loop = True

                else:
                    tmp_strings = []
                    for line in sorted(tmp):
                        tmp_strings += [item_sep.join(line)]
                    out_string += entry_start+entry_sep.join(tmp_strings)+entry_close

                    tmp = pointer[idx-1][0]
                    del pointer[idx]
    
        if fileformat == 'tex':
            out_string = out_string.replace('_', r'\_')
        tmpl = util.read_text_file(template) if template else '{0}'
        _write_file(keywords['filename'], tmpl.format(out_string), fileformat)

    def export(
            self,
            fileformat,
            sections=None,
            entries=None,
            entry_sep='',
            item_sep='',
            template='',
            **keywords
            ):
        """
        Export the wordlist to specific fileformats.

        Notes
        -----
        The difference between export and output is that the latter mostly
        serves for internal purposes and formats, while the former serves for
        publication of data, using specific, nested statements to create, for
        example, HTML or LaTeX files from the wordlist data.
        """

        self._export(
                fileformat,
                sections,
                entries,
                entry_sep,
                item_sep,
                template,
                **keywords
                )

    def tokenize(
            self,
            orthography_profile='',
            source="counterpart",
            target="tokens",
            column='graphemes',
            **keywords
            ):
        """
        Tokenize the data with help of orthography profiles.

        Parameters
        ----------
        ortho_profile : str (default='')
            Path to the orthographic profile used to convert and tokenize the 
            input data into IPA tokens.
        
        source : str (default="translation")
            The source data that shall be used for the tokenization procedures.
        
        target : str (default="tokens")
            The name of the target column that will be added to the wordlist.

        column : str (default="graphemes")
            Tokenization target.

        """
        self._tokenize(orthography_profile=orthography_profile, source=source, target=target,
                column=column, **keywords)

