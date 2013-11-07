# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-25 12:25
# modified : 2013-11-07 14:07
"""
Basic parser for text files in QLC format.
"""

__author__="Johann-Mattis List"
__date__="2013-11-07"

import os
import pickle
import codecs

from ..settings import rcParams
from ..read.qlc import read_qlc

from ..sequence.tokenizer import Tokenizer


class _QLCParser(object):
    """
    Basic class for the handling of text files in QLC format.

    """
    def __init__(
            self,
            filename,
            conf = '',
            ):
        
        # set the loaded var
        loaded = False

        # check for existing cache-directory
        if type(filename) == str:
            if os.path.isdir('__lingpy__'):
                # check for file extension
                if filename[-3:].lower() in ['qlc','csv']:
                    # split path in case it is over multiple dirs
                    fn = os.path.split(filename)[1]
                    path = os.path.join('__lingpy__',fn[:-3]+'bin')
                else:
                    path = ''
                if os.path.isfile(path):
                    if rcParams['verbose']: print("[i] Loading pickled object.")
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
        """
        Parse data regularly if the data has not been loaded from a pickled version.
        """

        # try to load the data
        internal_import = False

        # check whether it's a dictionary from which we load
        if type(filename) == dict:
            input_data = filename
            if 'filename' not in input_data:
                self.filename = rcParams['filename']
            internal_import = True
            
            # make check for correct input
            if len(input_data[0]) != len(input_data[1]):
                raise ValueError("[!] Wrong input format!")

        # check whether it's another wordlist-object
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
        
        # or whether the data is an actual file
        elif os.path.isfile(filename):
            input_data = read_qlc(filename)
            if filename[-3:].lower() in ['csv','qlc']:
                self.filename = filename[:-4]
            else:
                self.filename = filename
        
        # raise an error otherwise
        elif type(filename) == str:
            raise IOError("[ERROR] Input file '{0}' does not exist.".format(filename))
        else:
            raise TypeError("[ERROR] Unrecognized type for 'filename' arguemnt: {0}".format(type(filename).__name__))

        # load the configuration file
        if not conf:
            conf = os.path.join(rcParams['_path'],'data','conf','qlc.rc')

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
        
        # assign the data as attribute to the word list class. Note that we
        # need to check for the type here, but since numpy also offers integer
        # types, we don't check for type(x) == int, but instead use the
        # str.numeric-function that returns numeric values only if it is an
        # integer
        self._data = dict([(int(k),v) for k,v in input_data.items() if k != 0 and str(k).isnumeric()])

        # iterate over self._data and change the values according to the
        # functions (only needed when reading from file)
        if not internal_import:
            heads = sorted(self._header.items(),key=lambda x:x[1])
            for key in self._data:
                check = []
                for head,i in heads:
                    if i not in check:
                        self._data[key][i] = self._class[head](self._data[key][i])
                        check.append(i)

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
                    raise KeyError("[ERROR] The key {0} does not exist!".format(idx))

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

    def _pickle(self):
        """
        Store the current data in a pickled object.
        """
        if not os.path.isdir('__lingpy__'):
            os.mkdir('__lingpy__')
        
        # get the rest of the filename, important for cases in which a user
        # loads a wordlist inside a folder without having cded into it.
        fn = os.path.split(self.filename)[1]

        path = os.path.join('__lingpy__',fn+'.bin')
        out = open(path,'wb')
        d = {}
        for key,value in self.__dict__.items():
            if key not in  [
                    '_class',
                    ]:
                d[key] = value
        d['__date__'] = rcParams['timestamp']
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

    def _add_entries(
            self,
            entry,
            source,
            function,
            override=False,
            **keywords
            ):
        """
        Add new entry-types to the word list by modifying given ones.

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
            answer = input("[?] Datatype <{entry}> has already been produced, do you want to override? (y/n) ".format(entry=entry))
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
            # original data-storage of the wordlist
            elif hasattr(source,'keys') and hasattr(source,'values'): 
                
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

    def _tokenize(
            self,
            orthography_profile = '',
            source = "head",
            target = "tokens",
            conversion = 'graphemes',
            ** keywords
            ):
        """
        Tokenize the data with help of orthography profiles.
                
        """

        t = Tokenizer(orthography_profile)

        # else just return a Unicode grapheme clusters parse
        if target == 'tokens':
            function = lambda x: t.tokenize(x).split(' ')
        else:
            function = lambda x: t.tokenize(x)

        self.add_entries(
            target,
            source,
            function
            )

    def tokenize(
            self,
            orthography_profile = '',
            source = "counterpart",
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
        self._tokenize(orthography_profile, source, target, conversion,
            keywords)
