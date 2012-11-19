#! /usr/bin/env python3
"""
This module provides a basic class for the handling of word lists.

@date: 2012-11-08
@author: Johann-Mattis List
"""

import os
import builtins
import numpy as np

def _load_dict(infile):
    """
    Simple function only used to test the WordList class.
    """

    # read the data into a list
    data = []

    # open the file
    f = open(infile)

    for line in f:
        # ignore hashed lines
        if not line.startswith('#') and not line.startswith('@'):

            # mind to strip newlines
            data.append(line.strip('\n\r').split('\t'))
    
    # create the dictionary in which the data will be stored
    d = {}

    # check for first line, if a local ID is given in the header (or simply
    # "ID"), take this line as the ID, otherwise create it
    if data[0][0].lower() in ['local_id','localid']:
        local_id = True
    else:
        local_id = False

    # iterate over data and fill the dictionary (a bit inefficient, but enough
    # for the moment)
    i = 1
    for line in data[1:]:
        if local_id:
            d[int(line[0])] = line[1:]
        else:
            d[i] = line
            i += 1

    # assign the header to d[0]
    if local_id:
        d[0] = [x.lower() for x in data[0][1:]]
    else:
        d[0] = [x.lower() for x in data[0]]

    # return the stuff
    return d
            

class WordList(object):
    """
    Basic class for the handling of multilingual word lists.

    Parameters
    ----------
    data : dict
        A dictionary with consecutive integers as keys and lists as values.

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
            data,
            row = 'concept',
            col = 'doculect',
            conf = ''
            ):

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
        tmp = [line.strip('\n\r').split('\t') for line in open(conf)]
        
        # define two attributes, _alias, and _class which store the aliases and
        # the datatypes (classes) of the given entries
        self._alias,self._class = {},{}
        for name,cls,alias in tmp:
            # make sure the name itself is there
            self._alias[name.lower()] = name
            self._alias[name.upper()] = name

            self._class[name.lower()] = eval(cls)
            self._class[name.upper()] = eval(cls)

            # add the aliases
            for a in alias.split(','):
                self._alias[a.lower()] = name
                self._alias[a.upper()] = name
                self._class[a.lower()] = eval(cls)
                self._class[a.upper()] = eval(cls)

        # append the names in data[0] to self.conf to make sure that all data
        # is covered, even the types which are not specifically defined in the
        # conf file. the datatype defaults here to "str"
        for name in data[0]:
            if name.lower() not in self._alias:
                self._alias[name.lower()] = name.lower()
                self._class[name.lower()] = str
            if name.upper() not in self._alias:
                self._alias[name.upper()] = name.lower()
                self._alias[name.upper()] = str

        # add emtpy alias for empty strings
        self._alias[''] = ''

        # retrieve basic types for rows and columns from the word list
        try:
            rowIdx = [i for i in range(len(data[0])) if \
                    self._alias[data[0][i]] == row][0]
            colIdx = [i for i in range(len(data[0])) if \
                    self._alias[data[0][i]] == col][0]
        except:
            #print(
            #    "[!] Could not find row and col in configuration or input file!"
            #    )
            raise ValueError("[!] Could not find row and col in configuration or input file!")
            #print([self._alias[i] for i in data[0]],row,col)

        basic_rows = sorted(
                set(
                    [data[k][rowIdx] for k in data if k > 0]
                    )
                )
        basic_cols = sorted(
                set(
                    [data[k][colIdx] for k in data if k > 0]
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
        for key,value in [(k,v) for k,v in data.items() if k > 0]:
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
                    [self._alias[x] for x in data[0]],
                    range(len(data[0]))
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
        self._data = dict([(k,v) for k,v in data.items() if k != 0])

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
                        
    def __getitem__(self,idx):
        """
        Method allows quick access to the data by passing the integer key.
        """
        try:
            return self._cache[idx]
        except:
            pass

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
                pass

    def __getattr__(
            self,
            attr
            ):
        """
        Define how attributes are overloaded.
        """
        try:
            # get the right name
            attr = self._alias[attr]

            if attr == self._row_name:
                return self.rows
            elif attr == self._col_name:
                return self.cols
            elif attr in self._header:
                return self.getEntries(attr)
            else:
                raise AttributeError("%r object has no attribute %r" %
                        (type(self).__class__,attr))
        except:
             raise AttributeError("%r object has no attribute %r" %
                    (type(self).__class__,attr))

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


    def getDict(
            self,
            col = '',
            row = '',
            entry = '',
            ):
        """
        Function returns dictionaries of the cells matched by the indices.

        Parameters
        ----------
        idxA : string
            The first index evaluated by the method. It should reflect the name
            of either one of the rows or one of the columns.
        
        idxB : string (default = '')
            The second index evaluated by the method. It can be used to specify
            the datatype of the rows or columns selected.

        Notes
        -----
        Tobeadded

        Return
        ------
        entries : dict
            A dictionary of keys and values specifying the selected part of the
            data.
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
            print("[!] Neither rows nor columns are selected!")
       

    def getList(
            self,
            row='',
            col='',
            entry='',
            flat=False
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

        Notes
        -----
        Tobeadded

        Examples
        --------
        All examples make use of the small sample dataset

        Return
        ------
        data : list
            A list representing the selected part of the
            data.
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
            print("[!] Neither rows nor columns are selected!")
    
    def getEntries(
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

    def addEntries(
            self,
            entry,
            source,
            function,
            **args
            ):
        """
        Add new entry-types to the word list by modifying given ones.

        Parameters
        ----------
        entry : string
            A string specifying the name of the new entry-type to be added to the
            word list.

        source : string
            A string specifying the basic entry-type that shall be modified. 

        function : function
            A function which is used to convert the source into the target
            value.

        keywords : {dict list}
            A dictionary of keywords that are passed as parameters to the
            function.

        Notes
        -----
        This method can be used to add new entry-types to the data by
        converting given ones. 
        """

        # check whether the stuff is already there
        if entry in self._header:
            answer = input("[?] Datatype has already been produced, do you want to override?")
            if answer.lower() in ['y','yes']:
                pass
            else:
                pass
        else:
            # get the index of the source in self
            idx = self._header[source]            

            # get the new index into the header
            self._header[entry.lower()] = max(self._header.values())+1
            self._alias[entry.lower()] = entry.lower()
            self.header[entry.lower()] = self._header[entry.lower()]

            # iterate over the data and create the new entry
            for key in self:
                
                # get the source
                s = self[key][idx]

                # transform s
                t = function(s,**args)

                # add
                self[key].append(t)
        
    
    def getEtymDict(
            self,
            ref = "cogid",
            entry = ''
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
        ref = self._alias[ref]

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
            cogIdx = self._header[ref]
            
            # iterate over all data
            for key in self:
                cogid = self[key][cogIdx]
                colIdx = self.cols.index(self[key][self._colIdx])
                try:
                    self._etym_dict[ref][cogid]
                except:
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

    def getPaps(
            self,
            ref = 'cogid'
            ):
        """
        Function returns a list of present-absent-patterns of a given word list.

        Parameters
        ----------
        ref : string (default = "cogid")
            The reference entry which is used to store the cognate ids.
        """
        
        try:
            return self._cache['#paps#',ref]
        except:
            pass

        try:
            self._etym_dict[ref]
        except KeyError:
            print("[!] Could not find the specified reference.")
            return

        paps = {}
        # retrieve the values 
        for key,values in self._etym_dict[ref].items():
            paps[key] = []
            for value in values:
                if not value:
                    paps[key].append(1)
                else:
                    paps[key].append(0)
        
        self._cache['#paps#',ref] = paps
        
        return paps
