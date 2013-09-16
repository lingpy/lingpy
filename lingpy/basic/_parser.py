# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-25 12:25
# modified : 2013-09-16 16:58
"""
Basic parser for text files in QLC format.
"""

__author__="Johann-Mattis List"
__date__="2013-09-16"

import os
import pickle
import codecs

from ..settings import rcParams
from ..read.qlc import read_qlc

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
                    path = os.path.join('__lingpy__',filename[:-3]+'bin')
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
            self.filename = rcParams['filename']
            internal_import = True

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
        else:
            if not os.path.isfile(filename):
                raise IOError(
                        "[ERROR] Input file does not exist."
                        )
            else:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                lines = traceback.format_exception(exc_type, exc_value,
                    exc_traceback)
                raise ValueError('[ERROR] Could not parse the input file. {0}'.format(lines))

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
        
        # assign the data as attribute to the word list class
        self._data = dict([(k,v) for k,v in input_data.items() if k != 0 and type(k) == int])
        
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
        path = os.path.join('__lingpy__/',self.filename+'.bin')
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
