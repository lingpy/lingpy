#! /usr/bin/env python3
"""
This module provides specific classes which are used to handle linguistic
datasets.
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
    
    meta : dict
        A simple dictionary of key-value pairs storing some metadata of the
        original dataset.

    conf : string (default=None)
        A string defining the path to the configuration file. 

    Notes
    -----
    A word list is created from a data-dictionary and a meta-dictionary.
    The first contains the main data (keys and lists of values,
    corresponding to a column specified by key with id 0), and the latter
    contains additional data which was read in when parsing a qlc-formatted
    input file. If meta or conf are left undefined, the default values will
    be used.

    conf is the path to the configuration file. this file consists of 3
    columns, the first column is the internal name which is reserved for
    a column in the input data, the second column defines the datatype, and
    the third column defines aliases, separated by a comma

    When creating a wl instance, a specific data-structure is created
    consisting of 
    * self._data : the original data
    * self._dict : a dictionary representation of columns and cells
    * self._array : a flat representation of the data as an array
    * self._idx : an index which is needed to get the data of the array
    
    .. todo:: Add more documentation...
    """

    def __init__(
            self,
            data,
            meta = None,
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
                    )[0] + '/data/conf/wordlist.conf'

        # read the file defined by its path in conf
        tmp = [line.strip('\n\r').split('\t') for line in open(conf)]
        self.conf = {}
        for name,dt,tp,alias in tmp:
            for a in alias.split(','):
                self.conf[a] = [dt,tp,name]

                # the item which has "ROW" in the config file is used to align
                # all entries in a row
                if tp == 'ROW':
                    rows = name
                # the item which has "COL" as config type is used to align all
                # entries in a column
                if tp == 'COL':
                    cols = name

        # append the names in data[0] to self.conf to make sure that all data
        # is covered, even the types which are not specifically defined in the
        # conf file
        for name in data[0]:
            if name.lower() not in self.conf:
                self.conf[name.lower()] = ['str','opt',name.lower()]

        # retrieve basic types for rows and columns from the word list
        try:
            rowIdx = [i for i in range(len(data[0])) if \
                    self.conf[data[0][i].lower()][2] == rows][0]
            colIdx = [i for i in range(len(data[0])) if \
                    self.conf[data[0][i].lower()][2] == cols][0]
        except:
            print(
                "[!] Could not find ROW and COL in configuration or input file!"
                )

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
        self._data = data

        # the header stores the indices of the data in the original data
        # dictionary
        self._header = dict(
                zip(
                    [self.conf[x.lower()][2] for x in self._data[0]],
                    range(len(data[0]))
                    )
                )
                
                        
    def __getitem__(self,idx):
        """
        Method allows quick access to the data by passing the integer key.
        """
        try:
            # return full data entry as list
            return self._data[idx]
        except:
            try:
                # return data entry with specified key word
                return self._data[idx[0]][self._header[idx[1]]]
            except:
                pass

    def getDict(
            self,
            idxA,
            idxB = ''
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
        data : dict
            A dictionary of keys and values specifying the selected part of the
            data.
        """

        if not idxB:
            
            # if the index points to the rows
            if idxA in self.rows:
                return self._dict[idxA]

            # if it points to the columns
            if idxA in self.cols:
                data = {}

                for i,j in  [(self._data[i][self._rowIdx],i) for 
                        i in self._array[:,self.cols.index(idxA)] if i != 0]:
                    try:
                        data[i] += [j]
                    except:
                        data[i] = [j]
                return data

        else:

            if idxA in self.rows:
                data = {}
                for key,value in self._dict[idxA].items():
                    data[key] = []
                    for v in value:
                        data[key].append(self[v][self._header[idxB]])
                
                return data

            if idxA in self.cols:
                data = {}
                for i in self._array[:,self.cols.index(idxA)]:
                    
                    if i != 0:
                        # get the row name
                        m = self._data[i][self._rowIdx]

                        # get the relevant entry
                        c = self._data[i][self._header[idxB]]
                        try:
                            data[m] += [c]
                        except:
                            data[m] = [c]
                return data

    def getList(
            self,
            idxA,
            idxB='',
            flat=False
            ):
        """
        Function returns lists of the cells matched by the indices.

        Parameters
        ----------
        idxA : string
            The first index evaluated by the method. It should reflect the name
            of either one of the rows or one of the columns.
        
        idxB : string (default = '')
            The second index evaluated by the method. It can be used to specify
            the datatype of the rows or columns selected.
        flat : bool (default = False)
            Specify whether the returned list should be one- or
            two-dimensional, or whether it should contain gaps or not.

        Notes
        -----
        Tobeadded

        Return
        ------
        data : list
            A list specifying the selected part of the
            data.
        """

        if not idxB:
            
            # if the index points to the rows
            if idxA in self.rows:
                data = self._array[self._idx[idxA]]
                if flat:
                    return [i for i in data.flatten() if i != 0]
                else:
                    return data.tolist()

            # if it points to the columns
            if idxA in self.cols:
                data = self._array[:,self.cols.index(idxA)]

                if flat:
                    return [i for i in data if i != 0]
                else:
                    return data.tolist()

            # if the index points to specific entry classes itself
            if idxA in self._header:
                tmp_func = eval(self.conf[idxA][0])
                data = []
                for row in self._array:
                    tmp = [0 for i in row]
                    for i,cell in enumerate(row):
                        if cell != 0:

                            tmp[i] = tmp_func(self._data[cell][self._header[idxA]])
                        else:
                            if self.conf[idxA][0] == 'int':
                                tmp[i] = 0
                            else:
                                tmp[i] = ''

                    data.append(tmp)
                return data

        else:
            if idxA in self.rows:
                tmp = self._array[self._idx[idxA]]
                if flat:
                    tmp = [i for i in tmp.flatten() if i != 0]
                    data = [self._data[i][self._header[idxB]] for i in tmp]
                else:
                    data = tmp.tolist() 
                    for i,line in enumerate(tmp):
                        for j,cell in enumerate(line):
                            data[i][j] = self._data[cell][self._header[idxB]]
                return data

            elif idxA in self.cols:
                tmp = self._array[:,self.cols.index(idxA)]
                if flat:
                    tmp = [i for i in tmp.flatten() if i != 0]
                    data = [self._data[i][self._header[idxB]] for i in tmp]
                else:
                    data = []
                    for i in tmp:
                        if i != 0:
                            data.append(self._data[i][self._header[idxB]])
                        else:
                            data.append('-')

                return data

    def add(
            self,
            source,
            target,
            function,
            **keywords
            ):
        """
        Add new tables to the word list.

        Parameters
        ----------
        source : string
            A string specifying the basic values which shall be modified. 

        target : string
            A string specifiying the name of the new values to be added to the
            word list.

        function : function
            A function which is used to convert the source into the target
            value.

        keywords : dict
            A dictionary of keywords that are passed as parameters to the
            function.

        Notes
        -----
        tba
        """

        pass 


