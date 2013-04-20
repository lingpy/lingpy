# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-14 00:21
# modified : 2013-04-06 22:54
"""
This module provides a basic class for the handling of word lists.
"""

__author__="Johann-Mattis List"
__date__="2013-04-06"

import os
from datetime import date,datetime
import numpy as np
import pickle

# basic lingpy imports
from ..read.csv import read_qlc
from ..convert import *
from ..check.messages import FileWriteMessage
try:
    from ..algorithm.cython import cluster
    from ..algorithm.cython import misc
except:
    from ..algorithm.cython import _cluster as cluster
    from ..algorithm.cython import _misc as misc

# import ortho-parser
from ..sequence.orthography import OrthographyParser

class Wordlist(object):
    """
    Basic class for the handling of multilingual word lists.

    Parameters
    ----------
    filename : { string dict }
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
    
    .. todo:: Add more documentation...

    .. date:: 2012-11-08
    """

    def __init__(
            self,
            filename,
            row = 'concept',
            col = 'doculect',
            conf = ''
            ):

        # set the loaded var
        loaded = False

        # check for existing cache-directory
        if os.path.isdir('__lingpy__'):
            path = '__lingpy__/'+filename.replace('.csv','.bin')
            if os.path.isfile(path):
                # open the infile
                infile = open(path,'rb')

                # load the dictionary
                d = pickle.load(infile)

                # close the infile
                infile.close()

                # set the attributes
                for key,val in d.items():
                    setattr(self,key,val)
                
                # reset the class attribute
                self._class = dict([(key,eval(value)) for key,value in
                    self._class_string.items()])
                
                loaded = True
            else:
                loaded = False
        if not loaded:
            self._init_first(filename,row,col,conf)

    def _init_first(
            self,
            filename,
            row,
            col,
            conf
            ):

        # try to load the data
        try:
            input_data = read_qlc(filename)
            self.filename = filename.replace('.csv','')
        except:
            if type(filename) == dict:
                input_data = filename
                self.filename = 'lingpy-{0}'.format(str(date.today()))
            else:
                
            #if not input_data:
                raise ValueError('[i] Input data is not specified.')
            #else:
            #    input_data = filename
            #    self.filename = 'lingpy-{0}'.format(str(date.today()))

        # load the configuration file
        if not conf:
            conf = os.path.split(
                    os.path.dirname(
                        os.path.abspath(
                            __file__
                            )
                        )
                    )[0] + '/data/conf/wordlist.rc'

        # read the file defined by its path in conf
        tmp = [line.strip('\n\r').split('\t') for line in open(conf) if not
                line.startswith('#') and line.strip('\n\r')]
        
        # define two attributes, _alias, and _class which store the aliases and
        # the datatypes (classes) of the given entries
        self._alias,self._class,self._class_string,self._alias2 = {},{},{},{}
        for name,cls,alias in tmp:
            
            # make sure the name itself is there
            self._alias[name.lower()] = name
            self._alias[name.upper()] = name

            self._class[name.lower()] = eval(cls)
            self._class[name.upper()] = eval(cls)

            self._class_string[name.lower()] = cls
            self._class_string[name.upper()] = cls

            # add the aliases
            for a in alias.split(','):
                self._alias[a.lower()] = name
                self._alias[a.upper()] = name
                self._class[a.lower()] = eval(cls)
                self._class[a.upper()] = eval(cls)
                self._class_string[a.lower()] = cls
                self._class_string[a.upper()] = cls

            self._alias2[name] = sorted(set(alias.split(','))) + [name]

        # append the names in data[0] to self.conf to make sure that all data
        # is covered, even the types which are not specifically defined in the
        # conf file. the datatype defaults here to "str"
        for name in input_data[0]:
            if name.lower() not in self._alias:
                self._alias[name.lower()] = name.lower()
                self._class[name.lower()] = str
            if name.upper() not in self._alias:
                self._alias[name.upper()] = name.lower()
                self._class[name.upper()] = str

        # add emtpy alias for empty strings XXX why was that? I can't remember
        # why this was important XXX
        self._alias[''] = ''
        
        # retrieve basic types for rows and columns from the word list
        try:
            rowIdx = [i for i in range(len(input_data[0])) if \
                    self._alias[input_data[0][i]] == row][0]
            colIdx = [i for i in range(len(input_data[0])) if \
                    self._alias[input_data[0][i]] == col][0]
        except:
            raise ValueError("[!] Could not find row and col in configuration or input file!")

        basic_rows = sorted(
                set(
                    [input_data[k][rowIdx] for k in input_data if k != 0 and type(k) == int]
                    )
                )
        basic_cols = sorted(
                set(
                    [input_data[k][colIdx] for k in input_data if k != 0 and type(k) == int]
                    )
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
        for key,value in [(k,v) for k,v in input_data.items() if k != 0 and type(k) == int]:
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
        
        # the header stores the indices of the data in the original data
        # dictionary
        self.header = dict(
                zip(
                    [self._alias[x] for x in input_data[0]],
                    range(len(input_data[0]))
                    )
                )

        # now create a specific header which has all aliases
        self._header = dict([(k,v) for k,v in self.header.items()])
        
        # assign all aliases to the header
        for alias in self._alias:
            try:
                idx = self._header[self._alias[alias]]
                self._header[alias] = idx
            except:
                pass
        
        # assign the data as attribute to the word list class
        self._data = dict([(k,v) for k,v in input_data.items() if k != 0 and type(k) == int])

        # iterate over self._data and change the values according to the
        # functions
        heads = sorted(self._header.items(),key=lambda x:x[1])
        for key in self._data:
            check = []
            for head,i in heads:
                if i not in check:
                    self._data[key][i] = self._class[head](self._data[key][i])
                    check.append(i)

        # define a cache dictionary for stored data for quick access
        self._cache = {}

        # create entry attribute of the wordlist
        self.entries = sorted(set([b.lower() for a,b in self._alias.items() if b]))
                       
        # assign meta-data
        self._meta = {}
        for key in [k for k in input_data if type(k) != int]:
            self._meta[key] = input_data[key]

        # check for taxa in meta
        if 'taxa' not in self._meta:
            self._meta['taxa'] = self.taxa

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

    def __len__(self):
        """
        Length of a Wordlist is the number of counterparts.
        """
        return len(self._data)

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
                return self._meta[attr]
            except:
                raise AttributeError("%r object has no attribute %r" %
                        (self.__class__,attr))

    def __iter__(self):
        """
        Iteration is overloaded by iterating over all keys in the basic data.
        """

        return iter([key for key in self._data.keys()])

    def _clean_cache(self):

        """
        Function cleans the cache.
        """
        del self._cache
        self._cache = {}

    def _pickle(self):
        """
        Store the current data in a pickled object.
        """
        if not os.path.isdir('__lingpy__'):
            os.mkdir('__lingpy__')
        path = '__lingpy__/'+self.filename+'.bin'
        out = open(path,'wb')
        d = {}
        for key,value in self.__dict__.items():
            if key not in  [
                    '_class',
                    ]:
                d[key] = value
        d['__date__'] = str(datetime.today())
        pickle.dump(d,out)
        out.close()

    def pickle(self):
        """
        Store a dump of the data in a binary file.

        Notes
        -----
        The function creates a folder ``__lingpy__`` on your system containing a
        binary file called ``FILENAME.bin`` with ``FILENAME`` corresponding to
        the name of the original CSV-file. Instantiating the same
        :py:class:`~lingpy.basic.wordlist.Wordlist` instance again will first
        check for already compiled binary files and, if they are there, load
        them instead of the CSV-file.
        """

        self._pickle()

    def get_dict(
            self,
            col = '',
            row = '',
            entry = '',
            **keywords
            ):
        """
        Function returns dictionaries of the cells matched by the indices.

        Parameters
        ----------
        col : string (default='')
            The column index evaluated by the method. It should contain one of the
            values in the row of the
            :py:class:`~lingpy.basic.wordlist.Wordlist` instance, usually a
            taxon (language) name.
        
        row : string (default='')
            The row index evaluated by the method. It should contain one of the
            values in the row of the
            :py:class:`~lingpy.basic.wordlist.Wordlist` instance, usually a
            taxon (language) name.
        
        entry : string (default = '')
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
        The 'col' and 'row' keywords in the function are all aliased according
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
        the 'flat' keyword should be set to c{True}.
        
        Select all words (orthographical representation) for the concept
        "Harry"::
        
            >>> wl.get_list(concept="Harry",entry="words")
            [['hæri', 'haralt', 'gari', 'gari']]            
        
        Note that the values of the list that is returned are always
        two-dimensional
        lists, since it is possible that the original file contains synonyms
        (multiple words corresponding to the same concept). If one wants to
        have a flat representation of the entries, the 'flat' keyword should be
        set to c{True}::

            >>> wl.get_list(concept="Harry",entry="words",flat=True)
            ['hæri', 'haralt', 'gari', 'gari']
            
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

    def add_entries(
            self,
            entry,
            source,
            function,
            override = False,
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
        # check for emtpy entries etc.
        if not entry:
            print("[i] Entry was not properly specified!")
            return
        
        # check for override stuff, this causes otherwise an error message
        if entry not in self.header and override:
            return self.add_entries(entry,source,function,override=False)

        # check whether the stuff is already there
        if entry in self._header and not override:
            answer = input(
                    "[?] Datatype <{entry}> has already been produced, do "\
                    +"you want to override? ".format(
                        entry=entry
                        )
                    )
            if answer.lower() in ['y','yes']:
                keywords['override'] = True
                self.add_entries(entry,source,function,**keywords)
            else:
                print("[i] ...aborting...")
                return
        elif not override:

            # get the new index into the header
            # add a new alias if this is not specified
            if entry.lower() not in self._alias2:
                self._alias2[entry.lower()] = [entry.lower(),entry.upper()]
                self._alias[entry.lower()] = entry.lower()
                self._alias[entry.upper()] = entry.lower()

            # get the true value
            name = self._alias[entry.lower()]

            # get the new index
            newIdx = max(self._header.values()) + 1
            
            # change the aliassed header for each entry in alias2
            for a in self._alias2[name]:
                self._header[a] = newIdx

            self.header[name] = self._header[name]

            # modify the entries attribute
            self.entries = sorted(set(self.entries + [entry]))
            
            # check for multiple entries (separated by comma)
            if ',' in source:
                sources = source.split(',')
                idxs = [self._header[s] for s in sources]

                # iterate over the data and create the new entry
                for key in self:

                    # get the id line
                    s = self[key]

                    # transform according to the function
                    t = function(s,idxs)

                    # add the stuff to the dictionary
                    self[key].append(t)

            # if the source is a dictionary, this dictionary will be directly added to the
            # original data-storage of the wordlist
            elif type(source) == dict:
                
                for key in self:
                    s = source[key]
                    t = function(s)
                    self[key].append(t)
            
            else:
                # get the index of the source in self
                idx = self._header[source]            

                # iterate over the data and create the new entry
                for key in self:
                    
                    # get the source
                    s = self[key][idx]

                    # transform s
                    t = function(s,**keywords)

                    # add
                    self[key].append(t)
        
        elif override:

            # get the index that shall be replaced
            rIdx = self._header[entry.lower()]
            
            # check for multiple entries (separated by comma)
            if ',' in source:
                sources = source.split(',')
                idxs = [self._header[s] for s in sources]

                # iterate over the data and create the new entry
                for key in self:

                    # get the id line
                    s = self[key]

                    # transform according to the function
                    t = function(s,idxs)

                    # add the stuff to the dictionary
                    self[key][rIdx] = t

            # if the source is a dictionary, this dictionary will be directly added to the
            # original data-storage of the wordlist
            elif type(source) == dict:
                
                for key in self:
                    s = source[key]
                    t = function(s)
                    self[key][rIdx] = t

            else:
                # get the index of the source in self
                idx = self._header[source]            

                # iterate over the data and create the new entry
                for key in self:
                    
                    # get the source
                    s = self[key][idx]

                    # transform s
                    t = function(s,**keywords)

                    # add
                    self[key][rIdx] = t

        # clear the cache
        self._clean_cache()
    
    def get_etymdict(
            self,
            ref = "cogid",
            entry = '',
            loans = False
            ):
        """
        Return an etymological dictionary representation of the word list.

        Parameters
        ----------
        ref : string (default = "cogid")
            The reference entry which is used to store the cognate ids.

        entry : string (default = '')
            The entry-type which shall be selected.

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
            if not loans:
                for key in self:
                    cogid = self[key][cogIdx]
                    colIdx = self.cols.index(self[key][self._colIdx])

                    # assign new line if cogid is missing
                    if cogid not in self._etym_dict[ref]:
                        self._etym_dict[ref][cogid] = [0 for i in range(self.width)]
                    
                    # assign the values for the current session
                    try:
                        self._etym_dict[ref][cogid][colIdx] += [key]
                    except:
                        self._etym_dict[ref][cogid][colIdx] = [key]
            else:
                for key in self:
                    cogid = abs(self[key][cogIdx])
                    colIdx = self.cols.index(self[key][self._colIdx])

                    # assign new line if cogid is missing
                    if cogid not in self._etym_dict[ref]:
                        self._etym_dict[ref][cogid] = [0 for i in range(self.width)]
                    
                    # assign the values for the current session
                    try:
                        self._etym_dict[ref][cogid][colIdx] += [key]
                    except:
                        self._etym_dict[ref][cogid][colIdx] = [key]
        
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
            ref = 'cogid',
            entry = 'concept',
            missing = 0,
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
        
        etym_dict = self.get_etymdict(ref=ref,entry=entry)            

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
            taxa = 'taxa',
            concepts = 'concepts',
            cognates = 'cogid',
            threshold = 0.6,
            verbose = False,
            **keywords
            ):
        """
        Function calculates specific data.

        Parameters
        ----------
        data : str
            The type of data that shall be calculated.

        """
        # XXX take care of keywords XXX
        if data in ['distances','dst']:
            self._meta['distances'] = wl2dst(self,taxa,concepts,cognates)
        elif data in ['tre','nwk','tree']:
            if 'distances' not in self._meta:
                self.calculate('distances',taxa,concepts,cognates)
            if 'distances' not in keywords:
                keywords['distances'] = False
            if 'tree_calc' not in keywords:
                keywords['tree_calc'] = 'neighbor'

            self._meta['tree'] = matrix2tree(
                    self._meta['distances'],
                    self.taxa,
                    keywords['tree_calc'],
                    keywords['distances']
                    )

        elif data in ['groups','cluster']:
            if 'distances' not in self._meta:
                self.calculate('distances',taxa,concepts,cognates)
            self._meta['groups'] = matrix2groups(
                    threshold,
                    self.distances,
                    self.taxa
                    )
        else:
            return

        if verbose: print("[i] Successfully calculated {0}.".format(data))

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
                'ref'       : 'cogid',
                'entry'     : 'concept',
                'missing'   : 0,
                'filename'  : 'lingpy-{0}'.format(str(date.today())),
                'formatter' : 'concept',
                'tree_calc' : 'neighbor',
                'distances' : False,
                'ref'       : 'cogid',
                'threshold' : 0.6, # threshold for flat clustering
                'subset'    : False, # setup a subset of the data,
                'cols'      : False,
                'rows'      : False,
                'meta'      : self._meta,
                'entry'     : 'word',
                'taxa'      : False
                }
            
        # compare with keywords and add missing ones
        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        if fileformat in ['paps.nex','paps.csv']:
            paps = self.get_paps(
                    ref=keywords['ref'],
                    entry=keywords['entry'],
                    missing=keywords['missing']
                    )
            if fileformat == 'paps.nex':
                pap2nex(
                        self.cols,
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
        if fileformat == 'taxa':
            if hasattr(self,'taxa'):
                out = ''
                for taxon in self.taxa:
                    out += taxon + '\n'
                f = open(keywords['filename'] + '.taxa','w')
                f.write(out)
                f.close()
            else:
                raise ValueError(
                        "[i] Taxa are not available."
                        )
        
        # csv-output
        if fileformat == 'csv':
            
            # get the header line
            header = sorted(
                    [s for s in set(self._alias.values()) if s in self._header],
                    key = lambda x: self._header[x]
                    )
            header = [h.upper() for h in header]

            # get the data, in case a subset is chosen
            if not keywords['subset']:
                # write stuff to file
                wl2csv(
                        header,
                        self._data,
                        **keywords
                        )
            else:
                cols,rows = keywords['cols'],keywords['rows']

                if type(cols) not in [list,tuple,bool]:
                    raise ValueError(
                            "[i] Argument 'cols' should be list or tuple."
                            )
                if type(rows) not in [dict,bool]:
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

                # get the data
                out = {}

                for key,line in self._data.items():
                    if rows:
                        if eval(" and ".join(stmts)):
                            out[key] = [line[i] for i in indices]
                    else:
                        out[key] = [line[i] for i in indices]

                wl2csv(
                        header,
                        out,
                        **keywords
                        )
        
        # output dst-format (phylip)
        if fileformat == 'dst':

            # check for distances as keyword
            if 'distances' not in self._meta:
                self._meta['distances'] = wl2dst(self)
            
            # write data to file
            filename = keywords['filename']
            f = open(filename+'.'+fileformat,'w')
            out = matrix2dst(
                    self._meta['distances'],
                    self.taxa,
                    stamp=keywords['stamp']
                    )
            f.write(out)
            f.close()

            # display file-write-message
            FileWriteMessage(filename,fileformat).message('written')
        
        # output tre-format (newick)
        if fileformat in ['tre','nwk']: #,'cluster','groups']:
            
            # XXX bad line, just for convenience at the moment
            filename = keywords['filename']
            
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

            f = open(filename+'.'+fileformat,'w')
            f.write('{0}'.format(tree))
            f.close()

            FileWriteMessage(filename,fileformat).message('written')

        if fileformat in ['cluster','groups']:

            filename = keywords['filename']
            
            if 'distances' not in self._meta:
                
                self._meta['distances'] = wl2dst(self) # check for keywords

            if 'groups' not in self._meta:
                self._meta['groups'] = matrix2groups(
                        keywords['threshold'],
                        self._meta['distances'],
                        self.taxa
                        )

            f = open(filename+'.'+fileformat,'w')
            for taxon,group in sorted(self._meta['groups'].items(),key=lambda x:x[0]):
                f.write('{0}\t{1}\n'.format(taxon,group))
            f.close()

            FileWriteMessage(filename,fileformat).message('written')

        if fileformat in ['starling','star.csv']:

            # make lambda inline for data-check
            l = lambda x: ['-' if x == 0 else x][0]

            fileformat = 'starling_'+keywords['entry']+'.csv'

            f = open(keywords['filename']+'.'+fileformat,'w')
            if 'cognates' not in keywords:
                f.write('ID\tConcept\t'+'\t'.join(self.taxa)+'\n')
                for i,concept in enumerate(self.concepts):
                    lines = self.get_list(row=concept,entry=keywords['entry'])
                    for line in lines:
                        f.write(str(i+1)+'\t'+concept+'\t'+'\t'.join([l(t) for t in line])+'\n')
            else:
                f.write('ID\tConcept\t'+'\t'.join(['{0}\t COG'.format(t) for t
                    in self.taxa])+'\n')
                for i,concept in enumerate(self.concepts):
                    lines = self.get_list(row=concept,entry=keywords['entry'])
                    cogs = self.get_list(row=concept,entry=keywords['cognates'])
                    for j,line in enumerate(lines):
                        f.write(str(i+1)+'\t'+concept+'\t')
                        f.write('\t'.join('{0}\t{1}'.format(l(a),b) for a,b in
                            zip(line,cogs[j]))+'\n')
            f.close()
            FileWriteMessage(keywords['filename'],fileformat).message('written')

    def output(
            self,
            fileformat,
            **keywords
            ):
        """
        Write wordlist to file.

        Parameters
        ----------
        fileformat : {'csv', 'tre','nwk','dst', 'taxa', 'starling', 'paps.nex', 'paps.csv'}
            The format that is written to file. This corresponds to the file
            extension, thus 'csv' creates a file in csv-format, 'dst' creates
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

    def tokenize(
            self,
            ortho_profile = '',
            source = "counterpart",
            target = "tokens",
            ** keywords
            ):
        """
        Tokenize the data with help of orthography profiles.

        Parameters
        ----------
        ortho_profile : str (default='')
            Path to the orthographic profile used to convert and tokenize the
            input data into IPA tokens. If not specified, a simple Unicode
            grapheme parsing is carried out.

        source : str (default="counterpart")
            The source data that shall be used for the tokenization procedures.

        target : str (default="tokens")
            The name of the target column that will be added to the wordlist.
        
        Notes
        -----
        This is a shortcut to the extended
        :py:class:`~lingpy.basic.wordlist.Wordlist` class that loads data and
        automatically tokenizes it.
        """
        
        if os.path.exists(ortho_profile):
            ortho_path = ortho_profile
        else:
            ortho_path = os.path.split(
                    os.path.dirname(
                        os.path.abspath(
                            __file__
                            )
                        )
                    )[0] + '/data/orthography_profiles/' + ortho_profile
        
        # if the orthography profile does exist, carry out to tokenize the data
        if os.path.exists(ortho_path):
            op = OrthographyParser(ortho_path)

            # check for valid IPA parse
            if op.exists_multiple_columns():

                # tokenize the data, define specific output if target == 'tokens'
                # for direct lexstat input
                if target == 'tokens':
                    function = lambda x: op.graphemes_to_ipa(x).split(' ')[1:-1]
                else:
                    function = lambda x: op.graphemes_to_ipa(x)

                self.add_entries(
                        target,
                        source,
                        function
                        )
        
            
