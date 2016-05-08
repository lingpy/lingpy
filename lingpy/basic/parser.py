# *-* coding: utf-8 *-*
"""
Basic parser for text files in QLC format.
"""
from __future__ import unicode_literals, division, print_function
import os
import numpy as np
from collections import defaultdict

from six import text_type as str
from six import string_types

from lingpy.settings import rcParams
from lingpy.read.qlc import read_qlc
from lingpy import cache
from lingpy import util
from lingpy.util import confirm
from lingpy import log


class QLCParser(object):
    """
    Basic class for the handling of text files in QLC format.

    """
    @staticmethod
    def unpickle(filename):
        obj = cache.load(filename)
        obj._recreate_unpicklables()
        return obj

    def __init__(self, filename, conf=''):
        """
        Parse data regularly if the data has not been loaded from a pickled version.
        """
        self.log = log.get_logger()

        # try to load the data
        internal_import = False

        # check whether it's a dictionary from which we load
        if isinstance(filename, dict):
            input_data = filename
            if 'filename' not in input_data:
                self.filename = rcParams['filename']
            internal_import = True

            # make check for correct input, there was a bug with a wrong
            # evaluation which is hopefully fixed by now
            tmp_keys = [k for k in input_data if isinstance(k, int)]
            if len(input_data[0]) != len(input_data[tmp_keys[0]]):
                raise ValueError("[!] Wrong input format!")  # pragma: no cover
        # check whether it's another wordlist-object
        elif hasattr(filename, '_data') and hasattr(filename, '_meta'):
            input_data = dict(filename._data.items())
            input_data.update(filename._meta.items())
            input_data[0] = [a for a, b in sorted(
                filename.header.items(),
                key=lambda x:x[1],
                reverse=False)]
            internal_import = True
            self.filename = rcParams['filename']
        # or whether the data is an actual file
        elif isinstance(filename, string_types) and os.path.isfile(filename):
            input_data = read_qlc(filename)
            self.filename = filename
        # raise an error otherwise
        elif isinstance(filename, string_types):
            raise IOError("Input file '{0}' does not exist.".format(filename))
        else:
            raise TypeError("Unrecognized type for 'filename' argument: {0}".format(
                type(filename).__name__))

        # load the configuration file
        if not conf:
            conf = util.data_path('conf', 'qlc.rc')

        # read the file defined by its path in conf
        tmp = [line.split('\t') for line in util.read_config_file(conf)]

        # define two attributes, _alias, and _class which store the aliases and
        # the datatypes (classes) of the given entries
        self._alias, self._class, self._class_string, self._alias2 = {}, {}, {}, {}
        for name, cls, alias in tmp:
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

        # add empty alias for empty strings XXX why was that? I can't remember
        # why this was important XXX
        self._alias[''] = ''

        # the header stores the indices of the data in the original data
        # dictionary
        self.header = dict(
            zip([self._alias[x] for x in input_data[0]], range(len(input_data[0]))))

        # now create a specific header which has all aliases
        self._header = dict([(k, v) for k, v in self.header.items()])

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
        self._data = {
            int(k): v for k, v in input_data.items() if k != 0 and str(k).isnumeric()}
        # check for same length of all columns
        check_errors = ''
        for k, v in self._data.items():
            if len(v) != len(self.header):
                check_errors += 'Row {0} in your data contains {1} fields (expected {2})\n'.format(k, len(v), len(self.header))
        if check_errors:
            raise ValueError(check_errors)

        # iterate over self._data and change the values according to the
        # functions (only needed when reading from file)
        if not internal_import:
            heads = sorted(self._header.items(), key=lambda x: x[1])
            for key in self._data:
                check = []
                for head, i in heads:
                    if i not in check:
                        logstring = 'Problem with row {0} in col {1}, expected' +\
                            ' «{4}» as datatype but received «{3}» ' +\
                            ' (ROW: {2}, entry {5}).'
                        try:
                            self._data[key][i] = self._class[head](self._data[key][i])
                            check.append(i)
                        except KeyError: 
                            log.warn(
                                logstring.format(
                                    key,
                                    i,
                                    '|'.join([str(x) for x in self._data[key]]),
                                    self._data[key][i],
                                    self._class[head],
                                    head))
                        except ValueError:
                            log.warn(
                                logstring.format(
                                    key,
                                    i,
                                    '|'.join([str(x) for x in self._data[key]]),
                                    self._data[key][i],
                                    self._class[head],
                                    head))

        # create entry attribute of the wordlist
        self.entries = sorted(set([b.lower() for a, b in self._alias.items() if b]))

        # assign meta-data
        self._meta = {}
        for key in [k for k in input_data if type(k) != int]:
            self._meta[key] = input_data[key]

    def _recreate_unpicklables(self):
        """run `eval` on the string representations."""
        self.log = log.get_logger()
        self._class = {key: eval(value) for key, value in self._class_string.items()}

    def __getitem__(self, idx):
        """
        Method allows quick access to the data by passing the integer key.

        Parameters
        ----------
        idx : { int, str, tuple }
            The index which you pass to the method, which can be either an
            integer which will return the respective line of your parsed
            data, or a string, which allows to accept all values which are
            stored in the _meta-attribute of your parsed data, or a tuple,
            consisting of an integer key for the line and a string key for the
            respective field.

        Examples
        --------
        Load LingPy and the test_data function which gives us access to the
        files in the test dataset accompanying LingPy (we use a Wordlist
        object, but this works likewise with QLCParser, LexStat, and
        Alignments)::

            >>> from lingpy import *
            >>> from lingpy.tests.util import test_data
            >>> wl = Wordlist(test_data('KSL.qlc'))

        Get the first line in the data::
        
            >>> wl[1]
            ['Albanian', 'all', '1', 'gjithë', 'ɟiθ', ['ɟ', 'i', 'θ'], 4]

        Get the first line in the data and specify the value for the "doculect"
        column::
            >>> wl[1,'doculect']
            'Albanian'
        
        Get the attribute "doculect" which is stored in the _meta-attribute of
        the Wordlist object and contains all language names in the data::

            >>> wl.doculect
            ['Albanian', 'English', 'French', 'German', 'Hawaiian', 'Navajo', 'Turkish']
        
        Notes
        -----
        This method raises a KeyError if
        * you pass two indices and the first index is not a valid ID for any
            line in your data,
        * you pass one index and this index does neither correspond to a valid
            line ID in your data, nor to a valid key of the _meta-attribute of the
            Parser object.
        It returns None if you pass a valid index, but the column in your data
        does not exist.
        """
        if idx in self._data:
            return self._data[idx]
        if idx in self._meta:
            return self._meta[idx]
        if isinstance(idx, tuple) and len(idx) == 2:
            try:
                return self._data[idx[0]][self.header[self._alias[idx[1]]]]
            except KeyError:
                if idx[0] in self._data:
                    return
                else:
                    raise KeyError("No line with ID {0} specified could be found.".format(
                        idx[0])) 
        raise KeyError("No entry with the specified key {0} could be found".format(
            idx))


        #raise KeyError("[ERROR] The key {0} does not exist!".format(idx))

    def __len__(self):
        """
        Length of a Wordlist is the number of counterparts.
        """
        return len(self._data)

    def __iter__(self):
        """
        Iteration is overloaded by iterating over all keys in the basic data.
        """

        return iter([key for key in self._data.keys()])

    def pickle(self, filename=None):
        """
        Store the QLCParser instance in a pickle file.

        Notes
        -----
        The function stores a binary file called ``FILENAME.pkl`` with ``FILENAME``
        corresponding to the name of the original file in the
        `user cache dir <https://github.com/ActiveState/appdirs#some-example-output>`_
        for lingpy on your system.
        To restore the instance from the pickle call
        :py:meth:`~lingpy.basic.parser.QLCParser.unpickle`.
        """
        # we reset the _class attribute, because it may contain unpicklable stuff, like
        # `eval`ed lambdas.
        self._class = {}
        self.log = None
        cache.dump(self, filename or self.filename)
        # after pickling we have to recreate the attribute.
        self._recreate_unpicklables()

    def add_entries(
            self,
            entry,
            source,
            function,
            override=False,
            **keywords):
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
        self._add_entries(entry, source, function, override=override, **keywords)

    def _add_entries(
            self,
            entry,
            source,
            function,
            override=False,
            **keywords):
        # check for empty entries etc.
        if not entry:
            raise ValueError('Entry was not properly specified.')

        lentry = entry.lower()

        def _apply(key, s, *args, **kwargs):
            try:
                res = function(s, *args, **kwargs)
            except:
                raise ValueError('Could not convert item "{0}" (ID: {1}).'.format(s, key))
            if override:
                self._data[key][self._header[lentry]] = res
            else:
                self._data[key].append(res)

        # check for override stuff, this causes otherwise an error message
        if entry not in self.header and override:
            return self.add_entries(entry, source, function, override=False)

        # check whether the stuff is already there
        if entry in self._header and not override:
            if confirm(
                "Column <{entry}> already exists, do you want to override?".format(
                    entry=entry)):
                keywords['override'] = True
                return self.add_entries(entry, source, function, **keywords)
            return  # pragma: no cover

        if not override:
            # get the new index into the header
            # add a new alias if this is not specified
            if entry.lower() not in self._alias2:
                self._alias2[entry.lower()] = [entry.lower(), entry.upper()]
                self._alias[entry.lower()] = entry.lower()
                self._alias[entry.upper()] = entry.lower()

            # get the true value
            name = self._alias[entry.lower()]

            # get the new index
            newIdx = max(self._header.values()) + 1

            # change the aliased header for each entry in alias2
            for a in self._alias2[name]:
                self._header[a] = newIdx

            self.header[name] = self._header[name]

            # modify the entries attribute
            self.entries = sorted(set(self.entries + [entry.lower()]))

        # check for multiple entries (separated by comma)
        if ',' in source:
            sources = source.split(',')
            idxs = [self._header[s] for s in sources]

            # iterate over the data and create the new entry
            for key in self:
                _apply(key, self[key], idxs)
        # if the source is a dictionary, this dictionary will be directly added to the
        # original data-storage of the wordlist
        elif isinstance(source, dict):
            for key in self:
                _apply(key, source[key])
        else:
            # get the index of the source in self
            idx = self._header[source]
            for key in self:
                _apply(key, self[key][idx], **keywords)

class QLCParserWithRowsAndCols(QLCParser):
    def __init__(self, filename, row, col, conf):
        QLCParser.__init__(self, filename, conf=conf)

        try:
            self._row_name = self._alias[row]
            self._col_name = self._alias[col]
            rowIdx = self.header[self._alias[row]]
            colIdx = self.header[self._alias[col]]
        except KeyError:
            raise ValueError("Could not find row or col in configuration or input file!")

        def unique_sorted(idx, key):
            return sorted(
                set([self._data[k][idx] for k in self._data
                     if k != 0 and isinstance(k, int)]),
                key=key)

        # define rows and cols as attributes of the word list
        self.rows = unique_sorted(rowIdx, lambda x: ('%s' % x).lower())
        self.cols = unique_sorted(colIdx, lambda x: x.lower())

        # define height and width of the word list
        self.height = len(self.rows)
        self.width = len(self.cols)

        # row and column index point to the place where the data of the main
        # items is stored in the original dictionary
        self._rowIdx = rowIdx
        self._colIdx = colIdx

        # create a basic array which assigns ids for the entries in a starling manner.
        # first, find out, how many items (== synonyms) are there maximally for each row
        self._dict = defaultdict(lambda: defaultdict(list))
        for key, value in [(k, v) for k, v in self._data.items()
                           if k != 0 and str(k).isnumeric()]:
            self._dict[value[rowIdx]][value[colIdx]].append(key)

        # We must cast to a regular dict to make the attribute picklable.
        self._dict = dict(self._dict)

        # create the array by counting the maximal number of occurrences, store
        # the row names separately in a dictionary
        tmp_list = []
        self._idx = {}

        count = 0
        for k, d in self._dict.items():
            self._idx[k] = []

            # get maximal amount of "synonyms"
            for i in range(max([len(x) for x in d.values()])):
                tmp = []
                for j in range(self.width):
                    try:
                        tmp.append(d[self.cols[j]][i])
                    except:
                        tmp.append(0)
                self._idx[k] += [count]
                count += 1
                tmp_list += [tmp]

        self._array = np.array(tmp_list)

    def __getattr__(self, attr):
        """
        Define how attributes are overloaded.
        """
        if attr.startswith('_'):
            return self.__getattribute__(attr)

        if attr in self._alias:
            # get the right name
            nattr = self._alias[attr]

            if nattr == self._row_name:
                return self.rows

            if nattr == self._col_name:
                return self.cols

            if nattr in self._header:
                return self.get_entries(nattr)

        if attr in self._meta:
            return self._meta[attr]

        raise AttributeError("%r object has no attribute %r" % (self.__class__, attr))

    def get_entries(self, entry):
        """
        Return all entries matching the given entry-type as a two-dimensional list.

        Parameters
        ----------
        entry : string
            The entry-type of the data that shall be returned in tabular
            format.
        """
        if entry in self._header:
            entries = []
            for row in self._array:
                entries.append(
                    [self[cell][self._header[entry]] if cell != 0 else 0 for cell in row])
            return entries
