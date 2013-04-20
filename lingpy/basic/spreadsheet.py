"""
This module provides a basic class for reading in a simple spreadsheet (delimited text file) for concepts and words in a set of languages.
"""

__author__="Steven Moran"
__date__="2013-04"


import sys 
import unicodedata
from time import gmtime, strftime
from datetime import date,datetime

# don't use direct imports !!!
from ..sequence.ngram import *
from ..read.csv import *


class Spreadsheet:
    """
    Basic class for reading spreadsheet data that has been outputted into a deliminted format, e.g. tab.

    # workflow

    0. dump to delimited format
    1. csv2list(fileformat)
    2. pass this module arguments (header:linenumber, data:rownumber, default 0:ids, 1:concepts, 2-n:languages)
    2. col and row names (range or integer) - can tell the number of languages and concepts, etc.
    2. irrelevant order (loop over the dictionary and use alias dictionary; doculect / language / )
    2. header
    2. spreadsheet.rc **keywords, use aliases
    x. define for the output - keyword for output e.g. full rows, or rows >= n, also have to define what we have; want specific languages, specific cognate IDs
    x. black list parsing (no empty cells, etc.)
    x. separator for multientries for keywords in the output as "," as default, etc. (list of separators)
    x. parse as a list and do a type check
    x. then flip into harry potter format
    x. then flip into wordlist format
    x. then add tokenization / orthographic parsing

    # add stats to harry potter output

    """
    def __init__(self, 
                 filename,
                 fileformat = None,
                 dtype = None, 
                 comment = '#',
                 sep = '\t',
                 header = 0, # row of the header
                 concepts = 0, # column for the concepts
                 languages = [], # columns with language data -- TODO: must be ints
                 blacklist = "", # location of blacklist
                 conf = "" # spreadsheet .rc file
                 ):

        self.filename = filename
        self.fileformat = fileformat
        self.dtype = dtype
        self.comment = comment
        self.sep = sep
        self.header = header
        self.concepts = concepts
        self.languages = languages

        # create a 2D array and Unicode normalize its contents
        self.spreadsheet = csv2list(self.filename, self.fileformat, self.dtype, self.comment, self.sep)
        self._normalize()

        # given the header, concepts, and languages extract the data for processing
        # return a new matrix
        self.matrix = self._process_data()


    def _process_data(self):
        """
        Take spreadsheet input data and create matrix data given concept and language constraints.
        """
        
        # if no parameters, assume simple spreadsheet format
        if self.concepts == 0 and len(self.languages) == 0:
            print("[i] No concepts or languages specified. Assuming concepts in column 1 and all other cols language data.\n")
            return self.spreadsheet

        # else create a matrix given user-specified input
        matrix = []
        spreadsheet_header = self.spreadsheet[self.header] # get original header

        # create new header
        matrix_header = []
        matrix_header.append(spreadsheet_header[self.concepts])

        # TODO: this self.languages input must be int values
        for language in self.languages:
            matrix_header.append(spreadsheet_header[language])
        matrix.append(matrix_header)

        # append the concepts and words in languages and append the rows
        n = self.header+1
        for i in range(self.header+1, len(self.spreadsheet)):
            # check for concept; if missing skip row
            if self.spreadsheet[i][self.concepts] == "":
                print("[i] Missing concept in row "+str(i)+". Skipping the row.")
                continue

            # print(str(n), len(self.spreadsheet[i]))
            # n += 1

            row = []
            row.append(self.spreadsheet[i][self.concepts])
            for language in self.languages:
                # print("row", str(len(self.spreadsheet[i])), self.spreadsheet[i])                
                row.append(self.spreadsheet[i][language])
            matrix.append(row)

        return matrix


    def _normalize(self):
        """
        Function to Unicode normalize (NFD) cells in the matrix.
        """
        for i in range(0, len(self.spreadsheet)):
            for j in range(0, len(self.spreadsheet[i])):
                normalized_cell = unicodedata.normalize("NFD", self.spreadsheet[i][j])
                if not normalized_cell == self.spreadsheet[i][j]:
                    print("[i] Cell at <"+self.spreadsheet[i][j]+"> ["+str(i)+","+str(j)+"] not in Unicode NFD. Normalizing.")
                    self.spreadsheet[i][j] = normalized_cell


    def get_matrix_full_rows(self):
        """
        Create a 2D matrix from only the full rows in the spreadsheet.
        """
        full_row_matrix = []
        full_row_matrix.append(self.header)

        for row in self.matrix:
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
                    

    def print_matrix_stats(self):
        """
        Convenience function to get some stats data about the spreadsheet
        """
        total_entries = 0
        entries = []
        total_cells = len(self.matrix)*len(self.header)

        for header in self.header:
            entries.append(0)

        for row in self.matrix:
            for i in range(0, len(row)):
                if not row[i] == "":
                    total_entries += 1
                    entries[i] = entries[i] + 1
                
        print("total rows in matrix:", len(self.matrix))
        print("total cols in matrix:", len(self.header))
        print("total possible cells:", total_cells)
        print("total filled cells  :", str(total_entries), "("+str((total_entries*1.0)/total_cells*100)[:4]+"%)")
        print()
        print("total cells per column:")
        for i in range(0, len(self.header)):
            print(self.header[i], "\t", entries[i])

    
    def pprint(self):
        self.print_matrix(self.matrix)

    def print_matrix(self, matrix):
        """
        Print a matrix in tab delimited format
        """
        for i in range(0, len(matrix)):
            row = ""
            for j in range(0, len(matrix[i])):
                row += matrix[i][j]+"\t"
            row = row.rstrip("\t")
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

    def _output(self):
        filename = "lingpy-{0}".format(str(date.today()))
        output = open(filename, "w")
        for row in self.matrix:
            results = "\t".join(row)
            output.write(results+"\n")


    def output(self, **keywords):
        """
        Method to output the spreadsheet data into other formats.
        """

        """
        default:
        commas (",") == alternative transcriptions
        -- end of word comma == mistake
        backslash ("\") == two different words for the same thing
        periods (".") == identify affixes or units of meaning

        brackets ("[", "]") == 
        1a. after the word, e.g. mba [?] -- means "is it correct?"
        1b. after the word -- when an informant gives two different words
        2. within a word -- transcriber thinks segment is phonetic only
        3. after a word but not with "?" -- optional compounded word
        4. next to a word -- also optional compounded word

        # first do the data checks

        defaults = {
            "filename"  : "lingpy-{0}".format(str(date.today())),
            "subset"    : False,
            "blacklist" : None, # filename path
            "rows"      : 200, # minimum number of concepts
            "cols"      : 15 # minimum number of taxa
            }

        # add info to keywords
        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]
        """
        filename = "lingpy-{0}".format(str(date.today()))
        output = open(filename+".csv", "w")
        output.write("CONCEPT"+"\t"+"COUNTERPART"+"\t"+"COUNTERPART_DOCULECT"+"\n")
        for i in range(self.header+1, len(self.matrix)):
            for j in range(1, len(self.matrix[i])):
                output.write(self.matrix[i][0]+"\t"+self.matrix[i][j]+"\t"+self.matrix[self.header][j]+"\n")
            output.write("#\n")


if __name__=="__main__":
    # s = Spreadsheet("/Users/stiv/Dropbox/Projects/dogon/Moran_Dogon.comp.vocab.UNICODE.csv")
    # s = Spreadsheet("/Users/stiv/Dropbox/Projects/lingpy/test/spreadsheet_complex.tsv", concepts=2, languages=[3,5,6])
    # s = Spreadsheet("/Users/stiv/Dropbox/Projects/lingpy/test/spreadsheet_complex.tsv")
    # s = Spreadsheet("/Users/stiv/Dropbox/Projects/lingpy/test/Dogon_test.csv", concepts=7, languages=[17,18,19,20,22,23])


    s = Spreadsheet("/Users/stiv/Dropbox/Projects/lingpy/scripts/Dogon_100_DH.csv", concepts=2, languages=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22])
    s.output()

    # print(s.get_matrix_full_rows()) # bug with "0"

    # s.print_doculect_character_counts()

    # s.print_matrix_stats()
    # s.print_qlc_format()

    # frm = s.get_matrix_full_rows()
    # s.pprint(frm)
    # print(len(frm))
    # s.print_qlc_format(frm)

