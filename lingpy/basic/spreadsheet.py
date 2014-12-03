"""
Class for reading and tokenizing tabular data. e.g. csv data exported from a spreadsheet.
"""

__author__="Steven Moran"
__date__="2013-07"

# external imports
import sys 
import copy
import operator
import unicodedata
import collections
import operator
import regex as re
from time import gmtime, strftime

# internal imports
from ..sequence.tokenizer import *
from ..settings import rcParams
from ..sequence.ngram import *
from ..read.csv import *
from ..convert import *
from .. import util
from .. import log


class Spreadsheet:
    """
    Basic class to read tabular data (e.g. csv exported from spreadsheet), apply blacklist and pass to LingPy Wordlist.

    This module provides a basic class for reading in a simple spreadsheet (delimited text file) 
    for concepts and words in a set of languages. An optional blacklist can be used to remove / replace 
    characters, items, words, etc., from the input data. Some simple analyses of 

    Parameters
    ----------

    filename : str
        Name of the input file.

    fileformat : {None str}
        If not specified the file <filename> will be loaded. Otherwise, the
        fileformat is interpreted as the specific extension of the input file.

    dtype : {list}
        If not specified, all data will be loaded as strings. Otherwise, a
        list specifying the data for each line should be provided.

    comment : string (default="#")
        Comment character in the begin of a line forces this line to be
        ignored.

    sep : string (default = "\t")
        Specify the separator for the CSV-file.

    language_id : string (default = "NAME")
        A qualifier for columns in the input that contain language data.

    meanings : string (default = "CONCEPT")
        A qualifier for the concepts column in the input data.

    blacklist : string (default = "")
        Path and filename to the blacklist file.

    cellsep : string (default = ";")
        The character(s) indiciating the separtor between multiple linguistic forms in the same cell.

    skip_empty_cells : boolean (default = True)
        Flag for skipping empty cells in the input.


    """
    def __init__(self, 
                 filename,
                 fileformat = None, # required in lingpy/read/csv.py
                 dtype = None, # required in lingpy/read/csv.py
                 comment = '#', # comment marker
                 sep = '\t', # column separator
                 language_id = "NAME", # language column identifier
                 meanings = "CONCEPT", # explicit name of column containing concepts
                 blacklist = "", # filename path to blacklist file
                 cellsep = ';',
                 skip_empty_concepts = True,
                 **keywords
                 ):

        self.filename = filename
        self.fileformat = fileformat
        self.dtype = dtype
        self.comment = comment
        self.sep = sep
        self.language_id = language_id
        self.meanings = meanings
        self.blacklist = blacklist
        self.cellsep = cellsep
        self.skip_empty_concepts = skip_empty_concepts

        # set up matrix and normalize cell contents
        self.matrix = []
        self._init_matrix()
        self._normalize()

        # run blacklist replacements
        if self.blacklist: self._blacklist()

        # prepare dictionary for automatic pass-off to LingPy Wordlist
        self._prepare()


    def _init_matrix(self):
        """
        Create a 2D array from the CSV input and Unicode normalize its contents
        """
        # TODO: check if spreadsheet is empty and throw error
        spreadsheet = csv2list(
            self.filename, 
            self.fileformat, 
            self.dtype, 
            self.comment, 
            self.sep,
            strip_lines = False # this is of crucial importance, otherwise
            )

        # columns that have language data
        language_indices = []
        concept_id = 0

        # first row must be the header in the input; TODO: add more functionality
        header = spreadsheet[0] 

        log.info('%s' % header[0:10])

        for i, cell in enumerate(header):
            cell = cell.strip()
            log.info('%s' % cell)
            if cell == self.meanings:
                concept_id = i
            if self.language_id in cell:
                language_indices.append(i)

        matrix_header = []
        matrix_header.append(header[concept_id])        
        for i in language_indices:
            matrix_header.append(header[i].replace(self.language_id, "").strip())
        self.matrix.append(matrix_header)

        # append the concepts and words in languages and append the rows (skip header row)
        for i in range(1, len(spreadsheet)):
            matrix_row = []
            # if the concept cell is empty skip if flagged
            if spreadsheet[i][concept_id] == "" and self.skip_empty_concepts:
                continue
            for j in range(0, len(spreadsheet[i])):
                if j == concept_id or j in language_indices:
                    matrix_row.append(spreadsheet[i][j])
            self.matrix.append(matrix_row)

    def _normalize(self):
        """ 
        Function to Unicode normalize (NFD) cells in the input matrix.
        """
        for i in range(0, len(self.matrix)):
            for j in range(0, len(self.matrix[i])):
                normalized_cell = unicodedata.normalize("NFD", self.matrix[i][j])
                if not normalized_cell == self.matrix[i][j]:
                    log.debug("Cell at <"+self.matrix[i][j]+"> ["+str(i)+","+str(j)+"] not in Unicode NFD. Normalizing.")
                    self.matrix[i][j] = normalized_cell
    
    def _prepare(self, full_rows = False):
        """
        Prepare the spreadsheet for automatic pass-on to Wordlist.

        Notes
        -----
        We now assume that the matrix only contains concepts and counterparts.

        """

        # define a temporary matrix with full rows
        if not full_rows:
            matrix = self.matrix
        else:
            matrix = self.get_full_rows()

        # create the dictionary that stores all the data
        d = {}

        # iterate over the matrix
        idx = 1
        for i,line in enumerate(matrix[1:]):
            # only append lines that really work!
            if line:
                # get the concept
                concept = line[0].strip()
                if concept:
                    # get the rest
                    for j,cell in enumerate(line[1:]):
                        # get the language
                        language = matrix[0][j+1].replace(self.language_id,'').strip()
                        # get the counterparts
                        counterparts = [x.strip() for x in cell.split(self.cellsep)]
                        # append stuff to dictionary
                        for counterpart in counterparts:
                            if counterpart:
                                d[idx] = [concept,language,counterpart]
                                idx += 1

        # add the header to the dictionary
        d[0] = ["concept","doculect","counterpart"]

        # make the dictionary an attribute of spreadsheet
        self._data = dict([(k,v) for k,v in d.items() if k > 0])

        # make empty meta-attribute
        self._meta = dict(
                filename = self.filename
                )

        # make a simple header for wordlist import
        self.header = dict([(a,b) for a,b in zip(d[0],range(len(d[0])))])

    def _blacklist(self):
        """
        Remove anything in the spreadsheet that's specified in the blacklist file.

        """
        if not os.path.isfile(self.blacklist):
            log.warn("There is no blacklist specified at the follow file path location. Proceeding without blacklist.")
            return

        # loop through the blacklist file and compile the regexes
        rules = []
        replacements = []
        for line in util.read_config_file(self.blacklist, normalize='NFD'):
            rule, replacement = line.split(",,,") # black list regexes must be triple-comma delimited
            rule = rule.strip() # just in case there's trailing whitespace
            replacement = replacement.strip() # because there's probably trailing whitespace!
            rules.append(re.compile(rule))
            replacements.append(replacement)

        # blacklist the spreadsheet data - don't skip the header row, since
        # this may also contain blacklist information (as in Matthias' case)
        # doesn't apply to the header
        header = self.matrix[0]
        for i in range(1, len(self.matrix)):
            for j in range(0, len(self.matrix[i])):
                # skip removing anything from the concept column, e.g. ()'s
                if header[j] == self.meanings:
                    continue
                for k in range(0, len(rules)):
                    match = rules[k].search(self.matrix[i][j])
                    if not match == None:
                        match = re.sub(rules[k], replacements[k], self.matrix[i][j])
                        log.debug("Replacing ["+str(i)+","+str(j)+"] <"+self.matrix[i][j]+"> with <"+match+">.")
                        self.matrix[i][j] = match.strip()

    def pprint(self, matrix, delimit="\t"):
        """
        Convenice method to pretty print a 2D array.
        """
        for row in matrix:
            print("\t".join(str(x) for x in row))


    def analyze(self, *args):
        """ 
        This method analyses the character types, grapheme types and word types, 
        given a spreadsheet matrix.

        It returns a 2D matrix type of the type requested below.

        Notes
        -----
        *args: 

        "characters" = return a 2D matrix of chars (y) and counts in languages (x)
        "graphemes" = return a 2D matrix of graphemes (y) by counts in languages (x)
        "words" = return a 2D matrix of words by counts in languages

        """
        t = Tokenizer()

        char_types = collections.defaultdict(int)
        grapheme_types = collections.defaultdict(int)
        word_types = collections.defaultdict(int)

        chars_by_languages = collections.defaultdict(lambda : collections.defaultdict(int))
        graphemes_by_languages = collections.defaultdict(lambda : collections.defaultdict(int))
        words_by_languages = collections.defaultdict(lambda : collections.defaultdict(int))

        # deal with header and get language names
        header = self.matrix[0]

        # TODO put into English all counts across languages?
        # total_cells = len(self.matrix)*len(header)

        for i in range(1, len(self.matrix)):
            # skip empty rows
            if len(self.matrix[i]) == 0:
                continue
            # make sure rows aren't longer than the header row
            if len(self.matrix[i]) > len(header):
                log.error("You have a row (\#"+str(i)+") that is longer than your header. Exiting.")
                sys.exit(1)

            # process each cell for chars, graphemes, words and store the results
            for j in range(0, len(self.matrix[i])):
                cell = self.matrix[i][j].strip()
                if cell == "":
                    # TODO: integrate global verbosity
                    log.info("Missing cell")
                    continue
                # unicode characters
                for char in cell:
                    char_types[char] += 1
                    chars_by_languages[header[j]][char] += 1
                # unicode graphemes in qlc format
                grapheme_string = t.graphemes(cell)
                graphemes = grapheme_string.split()
                for grapheme in graphemes:
                    if not grapheme == "#":
                        grapheme_types[grapheme] += 1
                        graphemes_by_languages[header[j]][grapheme] += 1
                # word counts
                word_types[cell] += 1
                words_by_languages[header[j]][cell] += 1

        # sort the data types
        chars = sorted(char_types.items(), key=operator.itemgetter(0), reverse=False) 
        graphemes = sorted(grapheme_types.items(), key=operator.itemgetter(0), reverse=False) 
        words = sorted(word_types.items(), key=operator.itemgetter(0), reverse=False) 

        # deal with the user's request; if none, return tuple of 2D matrices
        results = []
        for arg in args:
            if arg == "characters":
                results.append(self._get_analysis("characters", header, chars, chars_by_languages))
            if arg == "graphemes":
                results.append(self._get_analysis("graphemes", header, graphemes, graphemes_by_languages))
            if arg == "words":
                results.append(self._get_analysis("words", header, words, words_by_languages))
        return results

    def _get_analysis(self, type, header, sorted_types, types_by_languages):
        """
        Private function to calculate and create a 2D matrix of the given types, e.g. a 
        languages by character types with x,y coordinates indicating character counts 
        per language.

        """
        # create and populate the matrices
        type_matrix = []
        type_header = [type]
        for item in header:
            if not item == self.meanings:
                type_header.append(item)
        type_matrix.append(type_header)

        # types are sort but the k,v are not in tuples
        for pair in sorted_types:
            result = []
            result.append(pair[0])
            for i in range(1, len(type_header)):
                result.append(types_by_languages[type_header[i]][pair[0]])
            type_matrix.append(result)
        return type_matrix

    def get_full_rows(self, matrix):
        """
        Create a 2D matrix from only the full rows in the spreadsheet.
        """
        full_row_matrix = []
        for row in matrix:
            is_full = 1
            for token in row:
                if token == "":
                    is_full = 0
            if is_full:
                full_row_matrix.append(row)
        return(full_row_matrix)

    def print_doculect_character_counts(self, doculects=1):
        for i in range(0, len(self.matrix)):
            print(self.matrix[i])
            for j in range(doculects, len(self.matrix[i])):
                if not self.matrix[i][j] == "":
                    print(self.matrix[i][j])
                    
    def stats(self):
        """
        Convenience function to get some stats data about the spreadsheet
        """
        total_entries = 0
        entries = []
        header = self.matrix[0]
        total_cells = len(self.matrix)*len(header)

        for item in header:
            entries.append(0)

        for row in self.matrix:
            for i in range(0, len(row)):
                if not row[i] == "":
                    total_entries += 1
                    entries[i] = entries[i] + 1
        print()
        print("### Simple matrix stats ###")
        print()
        print("# total cols in matrix:", len(header))
        print("# total rows in matrix:", len(self.matrix))
        print("# total possible cells:", total_cells)
        print("# total filled cells  :", str(total_entries), "("+str((total_entries*1.0)/total_cells*100)[:4]+"%)")
        print()
        print("# total cells per column:\n")
        for i in range(0, len(header)):
            print(str(entries[i]-1)+"\t"+str((entries[i]-1)/len(self.matrix)*100)+"\t"+header[i]) # do not include the header in count
        print()
    
    def pprint_matrix(self, delim="\t"):
        """
        Pretty print the spreadsheet matrix
        """
        for i in range(0, len(self.matrix)):
            row = ""
            for j in range(0, len(self.matrix[i])):
                row += self.matrix[i][j]+delim
            row = row.rstrip(delim)
            print(row)

    def print_qlc_format(self):
        """
        Print "simple" QLC format.
        """
        print("@input file: "+self.filename)
        print("@date: "+strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        print("#")
        print("LANGUAGE"+"\t"+"CONCEPT"+"\t"+"COUNTERPART")

        id = 0
        for i in range(1, len(self.matrix)):
            for j in range(1, len(self.matrix[i])):
                id += 1
                if self.matrix[i][j] == "":
                    row = str(id)+"\t"+self.header[j]+"\t"+self.matrix[i][0]+"\t"+"NaN"
                else:
                    row = str(id)+"\t"+self.header[j]+"\t"+self.matrix[i][0]+"\t"+self.matrix[i][j]
                print(row)        

    def _output(self, fileformat, **keywords):
        """
        Output the matrix into Harry Potter format.
        """
        defaults = dict(
                filename = "lingpy-{0}".format(_timestamp()),
                meta = self._meta
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]
        
        # use wl2csv to convert if fileformat is 'qlc'
        if fileformat in ['qlc','csv']:
            if fileformat == 'csv':
                log.deprecated('csv','qlc')
            wl2csv(
                    self.header,
                    self._data,
                    **keywords
                    )

    def output(
            self,
            fileformat,
            **keywords
            ):
        """
        Write Spreadsheet to different formats.
        """

        return self._output(fileformat,**keywords)
