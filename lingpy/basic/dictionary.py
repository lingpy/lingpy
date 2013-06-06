# author   : Peter Bouda
# email    : pbouda@cidles.eu
# created  : 2013-06-06 11:57
"""
This module provides a basic class for the handling of dictionaries.

TODO:
* What to do with dictionaries with more then 1 translation?

"""

__author__="Peter Bouda"
__date__="2013-06-06"

import os
import sys
from datetime import date,datetime
import numpy as np
import pickle
from operator import itemgetter
import pdb

# basic lingpy imports
from ..read.csv import read_qlc

class Dictionary():
    """
    Basic class for the handling of multilingual word lists.

    Parameters
    ----------
    filename : { string dict }
        The input file that contains the data. Otherwise a dictionary with
        consecutive integers as keys and lists as values with the key 0
        specifying the header.

    conf : string (default='')
        A string defining the path to the configuration file. 

    Notes
    -----
    A dictionary is created from a CSV file containing the data. Two keywords
    (head and translation) define, which of the dimensions of the original data
    should be used as heads and as translations of the dictionary content. A
    configuration file can be used to change basic names and aliases for the
    data being used, and the classes (data types) of the entries.

    """

    def __init__(
            self,
            filename,
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
            self._init_first(filename,conf)

    def _init_first(
            self,
            filename,
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
                raise ValueError('[i] Input data is not specified. Exception: {0}).'.format(sys.exc_info()))
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
                    )[0] + '/data/conf/dictionary.rc'

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

        # # define a cache dictionary for stored data for quick access
        # self._cache = {}

        # create entry attribute of the wordlist
        self.entries = sorted(set([b.lower() for a,b in self._alias.items() if b]))
                       
        # assign meta-data
        self._meta = {}
        for key in [k for k in input_data if type(k) != int]:
            self._meta[key] = input_data[key]


    def __getitem__(self,idx):
        """
        Method allows quick access to the data by passing the integer key.
        """
        try:
            # return full data entry as list
            out = self._data[idx]
            return out
        except:
            try:
                # return data entry with specified key word
                out = self._data[idx[0]][self._header[self._alias[idx[1]]]]
                return out
            except:
                try:
                    out = self._meta[idx]
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

            if nattr in self._header:
                return self.get_tuples([nattr])
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

    # def _clean_cache(self):

    #     """
    #     Function cleans the cache.
    #     """
    #     del self._cache
    #     self._cache = {}

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

    def get_tuples(self, columns = [ "head", "translation"]):
        """
        Return tuples from all entries for the given columns.

        Parameters
        ----------
        columns : list
            A list of the column names to extract from each entry.

        Returns
        -------
        entries : list
            A list of all the extracted columns. If there are more than
            one columns then each entry is a tuple. If there is only one column
            then each entry is a string.

        """

        # get the indices
        idxs = []
        for entry in columns:
            if entry in self._header:
                idxs.append(self._header[entry])

        entries = []
        if len(idxs) == 1:
            for row in self._data.values():
                entries.append(row[idxs[0]])
        elif len(idxs) > 1:
            mygetter = itemgetter(*idxs)
            for row in self._data.values():
                entries.append(tuple(mygetter(row)))

        return entries
