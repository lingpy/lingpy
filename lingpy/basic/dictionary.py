# author   : Peter Bouda
# email    : pbouda@cidles.eu
# created  : 2013-06-06 11:57
"""
This module provides a basic class for the handling of dictionaries.

"""

__author__="Peter Bouda"
__date__="2013-07-22"

import os
import sys
from datetime import date,datetime
import numpy as np
import pickle
import codecs
import re
from operator import itemgetter
import abc

# basic lingpy imports
from ..read.qlc import read_qlc
from ..settings import rcParams

# try:
#     from nltk.stem.snowball import SpanishStemmer
# except ImportError:
#    print(rcParams['W_missing_module'].format("nltk"))

# import tokenizer
from ..sequence.tokenizer import Tokenizer


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

        try:
            input_data = read_qlc(filename)
            if filename[-3:] in ['csv','qlc']:
                self.filename = filename[:-4]
            else:
                self.filename = filename
        except:
            if type(filename) == dict:
                input_data = filename
                self.filename = rcParams['filename'] 
            # if it's a wordlist object, add its basic parameters
            elif hasattr(filename,'_data') and hasattr(filename,'_meta'):
                input_data = dict(filename._data.items())
                input_data.update(filename._meta.items())
                input_data[0] = [a for a,b in sorted(
                    filename.header.items(),
                    key = lambda x:x[1],
                    reverse = False
                    )]
                internal_import = True
                self.filename = rcParams['filename'] 
            else:
                if not os.path.isfile(filename):
                    raise IOError(
                            "[i] Input file does not exist."
                            )
                else:
                    raise ValueError('[i] Could not parse the input file.')

        # load the configuration file
        if not conf:
            conf = os.path.join(rcParams['_path'],'data','conf','wordlist.rc')

        # read the file defined by its path in conf
        rcf = codecs.open(conf,'r','utf-8')
        tmp = [line.strip('\n\r').split('\t') for line in rcf if not
                line.startswith('#') and line.strip('\n\r')]
        rcf.close()
                
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

        # build a doculect->iso map for @doculect meta data
        self.doculect2iso = {}
        # do we have more than one doculect in the header?
        if type(self._meta["doculect"]) is list:
            for doculect in self._meta["doculect"]:
                doculect_entry = re.split(", ?", doculect)
                self.doculect2iso[doculect_entry[0]] = doculect_entry[1]
        else:
            doculect_entry = self._meta["doculect"]
            self.doculect2iso[doculect_entry[0]] = doculect_entry[1]

        # save ISO for heads and translations in list
        self.head_iso = []
        if type(self._meta["head_iso"]) is list:
            for iso in self._meta["head_iso"]:
                self.head_iso.append(iso)
        else:
            self.head_iso.append(self._meta["head_iso"])

        self.translation_iso = []
        if type(self._meta["translation_iso"]) is list:
            for iso in self._meta["translation_iso"]:
                self.translation_iso.append(iso)
        else:
            self.translation_iso.append(self._meta["translation_iso"])

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

    def add_entries(
            self,
            entry,
            source,
            function,
            override = False,
            **keywords
            ):
        """
        Add new entry-types to the dictionary by modifying given ones.

        Parameters
        ----------
        entry : string
            A string specifying the name of the new entry-type to be added to the
            dictionary.

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
            print(
                    "[?] Datatype <{entry}> has already been produced, ".format(entry=entry),
                    end = ''
                    )
            answer = input("do you want to override? (y/n) ")
            if answer.lower() in ['y','yes','j']:
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
            # original data-storage of the dictionary
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

    def tokenize(
            self,
            orthography_profile = '',
            source = "head",
            target = "tokens",
            conversion = 'graphemes',
            ** keywords
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

        conversion : str (default="graphemes")
            Tokenization target.


        Notes
        -----
        This is a shortcut to the extended
        :py:class:`~lingpy.basic.wordlist.Wordlist` class that loads data and
        automatically tokenizes it.
        
        """

        t = Tokenizer(orthography_profile)

        # else just return a Unicode grapheme clusters parse
        if target == 'tokens':
            function = lambda x: t.transform_rules(x).split(' ')
        else:
            function = lambda x: t.transform_rules(x)

        self.add_entries(
            target,
            source,
            function
            )
