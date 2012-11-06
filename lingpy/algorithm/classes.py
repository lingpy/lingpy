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
    if data[0][0].lower() in ['id','local_id','localid']:
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
        d[0] = data[0][1:]
    else:
        d[0] = data[0]

    # return the stuff
    return d

            

class WordList(object):
    """
    Basic class for the handling of multilingual word lists.

    A word list differs from a dictionary in so far as the keys of word lists
    are concepts...

    Todo
    ----

    Add more documentation...
    """

    def __init__(
            self,
            data,
            meta = None,
            conf = None
            ):
        """
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
        - self._data : the original data
        - self._dict : a dictionary representation of columns and cells
        - self._array : a flat representation of the data as an array
        - self._idx : an index which is needed to get the data of the array
        """

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
        Methods help to get quick access to the data in the word list.

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

        # carry out a type check of idx
        #t = type(idx)

        #if t == builtins.int:
        #    # return the data row if it's an integer key
        #    if idx in self._data:
        #        return self._data[idx]

        #elif t == builtins.str:
        #
        #    # return the dictionary if its a key to it
        #    if idx in self.rows:
        #        return self._dict[idx]

        #    ## return the data column if its a key to it
        #    #if idx in self.cols:
        #    #    data = [w for w in self._array[:,self.cols.index(idx)]:
        #    #        
        #    #    return zip(    
        #    #            self.rows,
        #    #            [w for w in self._array[:,self.cols.index(idx)] if w != 0]
        #    #            )

        #    # return a dictionary of the data if its a string-key
        #    try:
        #        if int(idx) in self._data:
        #            return dict(
        #                    zip(
        #                        [self.conf[x.lower()][2] for x in self._data[0]],
        #                        self._data[int(idx)]
        #                            )
        #                        )
        #    except:
        #        pass
        #
        ## if idx is a tuple, more specific return values can be defined. this
        ## is done by taking the first part of idx as the primary one and the
        ## second part as it's specifier
        #elif t == builtins.tuple:
        #    
        #    t2 = type(idx[0])

        #    if t2 == builtins.int:
        #        return self[idx[0]][self._header[idx[1]]]

        #    elif t2 == builtins.str:
        #        if idx[0] in self.rows:
        #            data = {}
        #            for key,value in self._dict[idx[0]].items():
        #                data[key] = [v for v in value]
        #                for i,v in enumerate(value):
        #                    data[key][i] = self._data[v][self._header[idx[1]]]

        #            return data

    def getDict(
            self,
            idxA,
            idxB = None
            ):
        """
        Function returns dictionaries of the cells matched by the indices.

        idxA should either occur in the rows of the word lists (e.g. be in the
        list of glosses), or in the columns (e.g. be in the list of taxa). idxB
        can be used to specify the values that shall returned, e.g. the entrys,
        the tokens, etc.
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
                for key,value in self._dict[idx[0]].items():
                    data[key] = []
                    for v in value:
                        data[key].append(self[v][self._header[idx[1]]])
                
                return data

            if idxA in self.cols:
                data = {}
                for i in self._array[:,self.cols.index(idx[0])]:
                    
                    if i != 0:
                        # get the row name
                        m = self._data[i][self._rowIdx]

                        # get the relevant entry
                        c = self._data[i][self._header[idx[1]]]
                        try:
                            data[m] += [c]
                        except:
                            data[m] = [c]
                return data

    def getList(
            self,
            idxA,
            idxB=None,
            flat=False
            ):
        """
        Return a list of the cells specified by maximally two parameters.
        """

        if not idxB:
            
            # if the index points to the rows
            if idxA in self.rows:
                out = self._array[self._idx[idxA]]
                if flat:
                    return [i for i in out.flatten() if i != 0]
                else:
                    return out.tolist()

            # if it points to the columns
            if idxA in self.cols:
                out = self._array[:,self.cols.index(idxA)]

                if flat:
                    return [i for i in out if i != 0]
                else:
                    return out.tolist()

            # if the index points to specific entry classes itself
            if idxA in self._header:
                tmp_func = eval(self.conf[idxA][0])
                out = []
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

                    out.append(tmp)
                return out

        else:

            if idxA in self.rows:
                tmp = self._array[self._idx[idxA]]
                if flat:
                    tmp = [i for i in tmp.flatten() if i != 0]
                    out = [self._data[i][self._header[idxB]] for i in tmp]
                else:
                    out = tmp.tolist() 
                    for i,line in enumerate(tmp):
                        for j,cell in enumerate(line):
                            out[i][j] = self._data[cell][self._header[idxB]]
                return out

            elif idxA in self.cols:
                tmp = self._array[:,self.cols.index(idxA)]
                if flat:
                    tmp = [i for i in tmp.flatten() if i != 0]
                    out = [self._data[i][self._header[idxB]] for i in tmp]
                else:
                    out = []
                    for i in tmp:
                        if i != 0:
                            out.append(self._data[i][self._header[idxB]])
                        else:
                            out.append('-')

                return out

