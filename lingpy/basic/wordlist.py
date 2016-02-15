# *-* coding: utf-8 *-* 
"""
This module provides a basic class for the handling of word lists.
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

import os
import numpy as np
import logging
from collections import defaultdict

from six import text_type as str

# basic lingpy imports
from ..convert.strings import matrix2dst, pap2nex, pap2csv, multistate2nex
from ..settings import rcParams
from .parser import QLCParserWithRowsAndCols
from .ops import wl2dst, wl2dict, renumber, clean_taxnames, calculate_data, \
        wl2qlc, tsv2triple, wl2multistate, coverage

from ..algorithm import clustering as cluster
from .. import util
from .. import log


def _write_file(filename, content, ext=None):
    if ext:
        filename = filename + '.' + ext
    util.write_text_file(filename, content)


class Wordlist(QLCParserWithRowsAndCols):
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
        QLCParserWithRowsAndCols.__init__(
            self, filename, row, col, conf or util.data_path('conf', 'wordlist.rc'))

        # setup other local temporary storage
        self._etym_dict = {}

        # check for taxa in meta
        if 'taxa' in self._alias:
            if self._alias['taxa'] not in self._meta:
                self._meta[self._alias['taxa']] = self.cols

    def __getitem__(self, idx):
        """
        Method allows quick access to the data by passing the integer key.
        """
        try:
            return self._get_cached(idx)
        except KeyError:
            try:
                # return data entry with specified key word
                self._cache[idx] = self._data[idx[0]][self._header[self._alias[idx[1]]]]
                return self._cache[idx]
            except:
                if idx in self._meta:
                    self._cache[idx] = self._meta[idx]
                    return self._cache[idx]

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
        assert not (row and col)

        if row or col:
            key = row or col
            attr = 'col' if col else 'row'

            if (key, entry) in self._cache:
                return self._cache[key, entry]

            if key not in getattr(self, attr + 's'):
                print("[!] The {0} you selected is not available!".format(attr))
                return

        if row:
            entries = self._dict[row]
            if entry:
                entries = {key: [self[i][self._header[entry]] for i in value]
                           for key, value in entries.items()}

            self._cache[row, entry] = entries
            return entries

        if col:
            data = defaultdict(list)
            for i, j in [(self[i][self._rowIdx], i)
                         for i in self._array[:, self.cols.index(col)] if i != 0]:
                data[i].append(j)

            entries = data
            if entry:
                entries = {key: [self[i][self._header[entry]] for i in value]
                           for key, value in entries.items()}

            self._cache[col, entry] = entries
            return entries

        for key in [k for k in keywords if k in self._alias]:
            if self._alias[key] == self._col_name:
                return self.get_dict(col=keywords[key], entry=entry)

            if self._alias[key] == self._row_name:
                return self.get_dict(row=keywords[key], entry=entry)

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

    def get_etymdict(
            self,
            ref="cogid",
            entry='',
            loans=False
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

        Returns
        -------
        etym_dict : dict
            An etymological dictionary representation of the data.

        Notes
        -----
        In contrast to the word-list representation of the data, an
        etymological dictionary representation sorts the counterparts according to
        the cognate sets of which they are reflexes. If more than one cognate
        ID are assigned to an entry, for example in cases of fuzzy cognate IDs
        or partial cognate IDs, the etymological dictionary will return one
        cognate set for each of the IDs.

        """
        
        # make an alias for the reference
        ref = (self._alias[ref],loans)

        # make converting function for loans
        if loans:
            f = lambda x: abs(x)
        else:
            f = lambda x: x

        # create an etymdict object
        if ref not in self._etym_dict:
            self._etym_dict[ref] = {}

            # get the index for the cognate id 
            cogIdx = self._header[ref[0]]
            
            # iterate over all data
            for key in self:
                cogids = self[key][cogIdx]
                colIdx = self.cols.index(self[key][self._colIdx])
                # check if data is not a list or tuple, if this is the case,
                # make it a fake-list, so we can treat it just as all the other
                # instances of fuzzy cognates (output is the same, though)
                if isinstance(cogids, (str, int, float)):
                    cogids = [cogids]
                    
                for cog in cogids:
                    cogid = f(cog)
                    # we initialize with zero here, since this corresponds to a
                    # missing entry in our data
                    if cogid not in self._etym_dict[ref]:
                        self._etym_dict[ref][cogid] = [0 for i in range(self.width)]
                    # assign new values for the current session
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
    
    def _output(self, fileformat, **keywords):
        """
        Internal function that eases its modification by daughter classes.
        """
        # check for stamp attribute
        keywords["stamp"] = getattr(self, '_stamp', '')

        # add the default parameters, they will be checked against the keywords
        util.setdefaults(
            keywords,
            cols=False,
            distances=False,
            entries=("concept", "counterpart"),
            entry='concept',
            fileformat=fileformat,
            filename=rcParams['filename'],
            formatter='concept',
            loans=False,
            meta=self._meta,
            missing=0,
            ref='cogid',
            rows=False,
            subset=False,  # setup a subset of the data,
            taxa='taxa',
            threshold=0.6,  # threshold for flat clustering
            tree_calc='neighbor')

        if fileformat in ['triple', 'triples', 'triples.tsv']:
            return tsv2triple(self, keywords['filename'] + '.' + fileformat)

        if fileformat in ['paps.nex', 'paps.csv']:
            paps = self.get_paps(
                ref=keywords['ref'], entry=keywords['entry'], missing=keywords['missing'])
            kw = dict(filename=keywords['filename'] + '.paps')
            if fileformat == 'paps.nex':
                kw['missing'] = keywords['missing']
                return pap2nex(self.cols, paps, **kw)
            return pap2csv(self.cols, paps, **kw)

        # simple printing of taxa
        if fileformat == 'taxa':
            assert hasattr(self, 'taxa')
            return util.write_text_file(keywords['filename'] + '.taxa', self.cols)

        # csv-output
        if fileformat in ['csv', 'qlc', 'tsv']:
            if fileformat == 'csv':
                log.deprecated('csv', 'qlc')

            # get the header line
            header = sorted(
                [s for s in set(self._alias.values()) if s in self._header],
                key=lambda x: self._header[x])
            header = [h.upper() for h in header]

            self._meta.setdefault('taxa', self.cols)

            # get the data, in case a subset is chosen
            if not keywords['subset']:
                # write stuff to file
                return wl2qlc(header, self._data, **keywords)

            cols, rows = keywords['cols'], keywords['rows']

            if not isinstance(cols, (list, tuple, bool)):
                raise ValueError("[i] Argument 'cols' should be list or tuple.")
            if not isinstance(rows, (dict, bool)):
                raise ValueError("[i] Argument 'rows' should be a dictionary.")

            # check for chosen header
            if cols:
                # get indices for header
                indices = [self._header[x] for x in cols]
                header = [c.upper() for c in cols]
            else:
                indices = [r for r in range(len(self.header))]

            if rows:
                stmts = []
                for key, value in rows.items():
                    if key == 'ID':
                        stmts += ["key " + value]
                    else:
                        idx = self._header[key]
                        stmts += ["line[{0}] ".format(idx) + value]

            log.debug("calculated what should be excluded")

            # get the data
            out = {}
            for key, line in self._data.items():
                if log.get_level() <= logging.DEBUG:
                    print(key)

                if rows:
                    if eval(" and ".join(stmts)):
                        out[key] = [line[i] for i in indices]
                else:
                    out[key] = [line[i] for i in indices]

            log.debug("passing data to wl2qlc")
            return wl2qlc(header, out, **keywords)
        
        # output dst-format (phylip)
        if fileformat == 'dst':
            # check for distances as keyword
            if 'distances' not in self._meta:
                self._meta['distances'] = wl2dst(self, **keywords)

            out = matrix2dst(self._meta['distances'], self.taxa, stamp=keywords['stamp'])
            return _write_file(keywords['filename'], out, fileformat)

        # output tre-format (newick)
        if fileformat in ['tre', 'nwk']:  # ,'cluster','groups']:
            if 'tree' not in self._meta:
                # check for distances
                if 'distances' not in self._meta:
                    self._meta['distances'] = wl2dst(self)
                # we look up a function to calculate a tree in the cluster module:
                tree = getattr(cluster, keywords['tree_calc'])(
                    self._meta['distances'], self.cols, distances=keywords['distances'])
            else:
                tree = self._meta['tree']

            return _write_file(keywords['filename'], '{0}'.format(tree), fileformat)

        if fileformat in ['cluster', 'groups']:
            if 'distances' not in self._meta:
                self._meta['distances'] = wl2dst(self)  # check for keywords

            if 'groups' not in self._meta:
                self._meta['groups'] = cluster.matrix2groups(
                    keywords['threshold'], self._meta['distances'], self.taxa)
            lines = []
            for taxon, group in sorted(self._meta['groups'].items(), key=lambda x: x[0]):
                lines.append('{0}\t{1}'.format(taxon, group))
            return _write_file(keywords['filename'], lines, fileformat)

        if fileformat in ['starling', 'star.csv']:
            # make lambda inline for data-check
            l = lambda x: ['-' if x == 0 else x][0]

            lines = []
            if 'cognates' not in keywords:
                lines.append('ID\tConcept\t' + '\t'.join(self.taxa))
                for i, concept in enumerate(self.concepts):
                    for line in self.get_list(row=concept, entry=keywords['entry']):
                        lines.append(
                            str(i + 1) + '\t' + concept + '\t' + '\t'.join(
                                [l(t) for t in line]))
            else:
                lines.append(
                    'ID\tConcept\t'+'\t'.join(['{0}\t COG'.format(t) for t in self.taxa]))
                for i, concept in enumerate(self.concepts):
                    cogs = self.get_list(row=concept, entry=keywords['cognates'])
                    for j, line in enumerate(self.get_list(row=concept, entry=keywords['entry'])):
                        part = '\t'.join('{0}\t{1}'.format(l(a), b) for a, b in zip(line, cogs[j]))
                        lines.append(str(i + 1) + '\t' + concept + '\t' + part)

            return _write_file(
                keywords['filename'], lines, 'starling_' + keywords['entry'] + '.csv')

        if fileformat == 'multistate.nex':
            if not keywords['filename'].endswith('.multistate.nex'):
                keywords['filename'] += '.multistate.nex'

            matrix = wl2multistate(self, keywords['ref'])
            return multistate2nex(self.taxa, matrix, keywords['filename'])

        if fileformat == 'separated':
            if not os.path.isdir(keywords['filename']):
                os.mkdir(keywords['filename'])

            for l in self.cols:
                lines = [''] if 'ignore_keys' in keywords else ['ID\t']
                lines[0] += '\t'.join([x.upper() for x in keywords['entries']])
                for key in self.get_list(col=l, flat=True):
                    line = [] if 'ignore_keys' in keywords else [key]
                    for entry in keywords['entries']:
                        tmp = self[key, entry]
                        if isinstance(tmp, list):
                            tmp = ' '.join([str(x) for x in tmp])
                        line += [tmp]
                    lines.append('\t'.join('{0}'.format(x) for x in line))
                _write_file('{0}/{1}'.format(keywords['filename'], l), lines, 'tsv')

    def output(self, fileformat, **keywords):
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
        return self._output(fileformat, **keywords)
    
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
        if not sections:
            if fileformat == 'txt':
                sections = dict(
                    h1=('concept', '\n# Concept: {0}\n'),
                    h2=('cogid', '## Cognate-ID: {0}\n'))
            elif fileformat == 'tex':
                sections = dict(
                    h1=('concept', r'\section{{Concept: ``{0}"}}' + '\n'),
                    h2=('cogid', r'\subsection{{Cognate Set: ``{0}"}}' + '\n'))
            elif fileformat == 'html':
                sections = dict(
                    h1=('concept', '<h1>Concept: {0}</h1>'),
                    h2=('cogid', '<h2>Cognate Set: {0}</h2>'))

        if not entries:
            if fileformat == 'txt':
                entries = [('language', '{0} '), ('ipa', '{0}\n')]
            elif fileformat == 'tex':
                entries = [('language', '{0} '), ('ipa', '[{0}]' + '\n')]
            elif fileformat == 'html':
                entries = [('language', '{0}&nbsp;'), ('ipa', '[{0}]\n')]
        
        util.setdefaults(keywords, filename=rcParams['filename'])

        # get the temporary dictionary
        out = wl2dict(self, sections, entries, exclude)
    
        # assign the output string
        out_string = ''
    
        # iterate over the dictionary and start to fill the string
        for key in sorted(out, key=lambda x: str(x).lower()):
            # write key to file
            out_string += key[1]
            
            # reassign tmp
            tmp = out[key]

            # set the pointer and the index
            pointer = {0: [tmp, sorted(tmp.keys())]}

            while True:
                idx = max(pointer.keys())

                # check for type of current point
                if isinstance(tmp, dict):
                    if pointer[idx][1]:
                        next_key = pointer[idx][1].pop()
                        out_string += next_key[1]
                        tmp = pointer[idx][0][next_key]
                        if isinstance(tmp, dict):
                            pointer[idx + 1] = [tmp, sorted(tmp.keys())]
                        else:
                            pointer[idx + 1] = [tmp, tmp]
                    else:
                        del pointer[idx]
                        if idx == 0:
                            break
                else:
                    tmp_strings = []
                    for line in sorted(tmp):
                        tmp_strings += [item_sep.join(line)]
                    out_string += entry_start + entry_sep.join(tmp_strings) + entry_close
                    tmp = pointer[idx - 1][0]
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
            **keywords)

    def coverage(self, stats='absolute'):
        """
        Function determines the coverage of a wordlist.
        """
        cov = coverage(self)

        if stats == 'absolute':
            return cov
        if stats == 'ratio':
            return dict([(a, b / self.height) for a, b in cov.items()])
        if stats == 'mean':
            return sum([a / self.height for a in cov.values()]) / self.width
    

