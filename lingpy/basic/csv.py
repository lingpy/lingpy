"""
Class for reading and tokenizing tabular data exported from a spreadsheet.
"""

__author__="Steven Moran"
__date__="2013-07"

# external imports
import sys 
import copy
import codecs
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


class Csv:
    """
    Basic class for reading tabular data exported from a spreadsheet.

    This module provides a basic class for reading in a simple spreadsheet (delimited text file) for concepts and words in a set of languages.

    Parameters
    ----------

    Notes
    -----
    TODO: how about a "config file" variable where the parameters below can be specified?

    """
    def __init__(self, 
                 filename,
                 fileformat = None, # req'd in read.csv
                 dtype = None, # req' in read.csv 
                 comment = '#', # comment marker
                 sep = '\t', # column separator
                 language_id = "NAME", # language column identifier
                 # lang_cols = [], # int specified columns for languages in a spreadsheet (index starts at 0)
                 meanings = "CONCEPT", # explicit name of column containing concepts
                 blacklist = "", # filename path to blacklist file
                 conf = "", # spreadsheet .rc file
                 cellsep = ';', # cell separator, separates forms in the same cell; delimited with pipe "|"
                 verbose = False,
                 skip_empty_concepts = True,
                 profiles = {}, # file path and names to orthography profiles
                 **keywords
                 ):

        self.filename = filename
        self.fileformat = fileformat
        self.dtype = dtype
        self.comment = comment
        self.sep = sep
        self.language_id = language_id
        # self.lang_cols = lang_cols
        self.meanings = meanings
        self.blacklist = blacklist
        self.conf = conf
        self.cellsep = cellsep
        self.skip_empty_concepts = skip_empty_concepts
        self.profiles = profiles
        self.verbose = verbose

        # set up matrix
        self.matrix = []
        self._init_matrix()
        self._normalize()

        # if self.blacklist: self._blacklist()
        if self.profiles: 
            self.tokenized_matrix = self._tokenize(profiles)

        # self._prepare()
        self._prepare_tokenized()

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

        if rcParams['verbose']: print(header[0:10])
        
        for i, cell in enumerate(header):
            cell = cell.strip()
            if rcParams['verbose']: print(cell)
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

        """
        matrix_header = []
        matrix_header.append(header[self.concepts])        
        for i in language_indices:
            matrix_header.append(header[i].replace(self.language_id, "").strip())
        self.matrix.append(matrix_header)

        # append the concepts and words in languages and append the rows
        for i in range(1, len(spreadsheet)): # skip the header row
            matrix_row = [] # collect concepts and languages to add to matrix
            temp = []
            for j in range(0, len(spreadsheet[i])):
                if j == self.concepts:
                    matrix_row.append(spreadsheet[i][j])
                if j in language_indices:
                    temp.append(spreadsheet[i][j])
            for item in temp:
                print(item)
                matrix_row.append(item)
            self.matrix.append(matrix_row)
            """

    def _normalize(self):
        """ 
        Function to Unicode normalize (NFD) cells in the matrix.
        """
        for i in range(0, len(self.matrix)):
            for j in range(0, len(self.matrix[i])):
                normalized_cell = unicodedata.normalize("NFD", self.matrix[i][j])
                if not normalized_cell == self.matrix[i][j]:
                    if self.verbose:
                        print("[i] Cell at <"+self.matrix[i][j]+"> ["+str(i)+","+str(j)+"] not in Unicode NFD. Normalizing.")
                    self.matrix[i][j] = normalized_cell
    
    def _prepare(self,full_rows = False):
        """
        Prepare the spreadsheet for automatic pass-on to Wordlist.
        """
        # XXX we now assume that the matrix is 'normalized',i.e. that it only
        # contains concepts and counterparts, in later versions, we should make
        # this more flexible by adding, for example, also proto-forms, or
        # cognate ids

        # define a temporary matrix with full rows
        if not full_rows:
            # matrix = self.matrix
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

        Notes
        -----

        TODO: make parameter **kwargs
        """
        if not os.path.isfile(self.blacklist):
            if rcParams['verbose']:
                print("[i] There is no blacklist specified at the follow file path location. Proceeding without blacklist.")
            return

        blacklist_file = codecs.open(self.blacklist, "r",'utf-8')
        # loop through the blacklist file and compile the regexes
        rules = []
        replacements = []
        for line in blacklist_file:
            line = line.strip()
            # skip any comments
            if line.startswith("#") or line == "":
                continue
            line = unicodedata.normalize("NFD", line)
            rule, replacement = line.split(",,,") # black list regexes must be triple-comma delimited
            rule = rule.strip() # just in case there's trailing whitespace
            replacement = replacement.strip() # because there's probably trailing whitespace!
            rules.append(re.compile(rule))
            replacements.append(replacement)
        blacklist_file.close()

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
                        if self.verbose:
                            print("[i] Replacing <"+self.matrix[i][j]+"> ["+str(i)+","+str(j)+"] with <"+match+">.")
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
                print("[i] You have a row (\#"+str(i)+") that is longer than your header. Exiting.")
                sys.exit(1)

            # process each cell for chars, graphemes, words and store the results
            for j in range(0, len(self.matrix[i])):
                cell = self.matrix[i][j].strip()
                if cell == "":
                    # TODO: integrate global verbosity
                    # if self.verbose:
                        # print("[i] Missing cell")
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


    def _tokenize(self, profiles):
        """ 
        Takes a dictionary or list of kwargs that specifies "column name" and "path to orthography profile".
        This function lets a user specify different orthographic profiles per column in a spreadsheet that's 
        already been semantically aligned.
        """
        tokenized_matrix = copy.deepcopy(self.matrix)

        # get orthography parsers and store them in dict
        ops = {}
        for k, v in profiles.items():
            ops[k] = Tokenizer(v) #, debug=True)

        # TODO: identify the column and orthographically parse it in the matrix

        # apply tokenization
        header = tokenized_matrix[0]
        for i in range(1, len(tokenized_matrix)):
            for j in range(0, len(tokenized_matrix[i])):
                cell = tokenized_matrix[i][j].strip()
                if cell == "":
                    continue

                # split up tokens by the cell separtor and grapheme-ically parse them
                tokens = cell.split(self.cellsep)
                tokenized_tokens = []
                for token in tokens:
                    token = token.strip()
                    if header[j] in ops:
                        tokenized_form = ops[header[j]].transform_rules(token).strip() # .replace("#", "").strip()
                        if tokenized_form == "":
                            tokenized_tokens.append("FAIL")
                        else:
                            tokenized_tokens.append(tokenized_form)
                            tokenized_matrix[i][j] = " \\\\ ".join(tokenized_tokens)

                    # original
                    # tokenized_matrix[i][j] = ops[header[j]].parse_graphemes(cell).strip() # .replace("#", "").strip()
                    # print(ops[header[j]].parse_graphemes(cell).strip()) # .replace("#", "").strip()
                    else:
                        if header[j] in ops:
                            tokenized_form = ops[header[j]].transform_rules(cell).strip() # .replace("#", "").strip()
                            if tokenized_form == "":
                                tokenized_matrix[i][j] = "FAIL"
                            else:
                                tokenized_matrix[i][j] = tokenized_form                    
        return tokenized_matrix

    def _prepare_tokenized(self, full_rows=False):
        """
        # Prepare the spreadsheet(s) for automatic pass-on to Wordlist.
        # We now assume that the matrix contains only concepts and counterparts, 
        # TODO: in later versions, we should make this more flexible by adding, 
        # for example, also proto-forms, or cognate ids, etc.
        """

        """
        # define temp matrices with full rows if specified
        if not full_rows:
            matrix = self.matrix
            tokenized_matrix = self.tokenized_matrix
        else:
            matrix = self.get_full_rows(self.matrix)
            tokenized_matrix = self.get_full_rows(self.tokenized_matrix)
            """
        matrix = self.matrix
        tokenized_matrix = self.tokenized_matrix

        # store and process the data
        wordlist = {}
        id = 1
        for i in range(1, len(tokenized_matrix)):
            # get concept
            if tokenized_matrix[i]:
                concept = tokenized_matrix[i][0].strip()
                # proceed only if there is a concept
                if concept:
                    for j in range(1, len(tokenized_matrix[i])):
                        language = tokenized_matrix[0][j].replace(self.language_id, "").strip()
                        counterparts = [x.strip() for x in tokenized_matrix[i][j].split(self.cellsep)]
                        # if a tokenized matrix, gather the tokenized forms
                        tokenized_counterparts = []
                        if self.tokenized_matrix:
                            tokenized_counterparts = [y.strip() for y in tokenized_matrix[i][j].split(self.cellsep)]

                        if not len(counterparts) == len(tokenized_counterparts):
                            print(concept)
                            print(self.matrix[i])
                            print(self.tokenized_matrix[i])
                            print(self.matrix[i][j])
                            print(len(counterparts), counterparts)
                            print(self.tokenized_matrix[i][j])
                            print(len(tokenized_counterparts), tokenized_counterparts)
                            print()


                        # append stuff to dictionary
                        for k in range (0, len(counterparts)):
                            if self.tokenized_matrix:
                                if counterparts[k]:
                                    # wordlist[id] = [concept, language, counterparts[k], tokenized_counterparts[k]]
                                    wordlist[id] = [concept, language, counterparts[k]]
                                    id += 1
                            else:
                                if counterparts[k]:
                                    wordlist[id] = [concept, language, counterparts[k]]
                                    id += 1

        # add the header to the wordlist dictionaries
        if self.tokenized_matrix:
            wordlist[0] = ["concept", "doculect", "counterpart", "tokens"]
        else:
            wordlist[0] = ["concept", "doculect", "counterpart"]

        # make the dictionary an attribute of this spreadsheet object
        self._data = dict([(k, v) for k, v in wordlist.items() if k > 0])

        # make empty meta-attribute
        self._meta = dict(
                filename = self.filename
                )

        # make a simple header for wordlist import
        self.header = dict([(a, b) for a, b in zip(wordlist[0], range(len(wordlist[0])))])


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
                print(rcParams['W_deprecation'].format('csv','qlc'))
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


