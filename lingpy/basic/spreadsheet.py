# -*- coding: utf-8 -*-
"""
This module provides a basic class for reading in a simple spreadsheet 
(delimited text file) for concepts and words in a set of languages.
"""

import csv
import sys 
import unicodedata
from time import gmtime, strftime

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

"""


class Spreadsheet:
    """
    Basic class for reading CSV files and outputting QLC format.

    Parameters
    ----------
    filename : str
       Path to the file to be loaded.

    delimiter : str
       The delimiter used in the CSV file. Default it tab.

    output_dir : str
       The output directory to write data to.

    Notes
    -----
    Here will be some notes.

    .. todo:: finish documentation; what to do about cell encapsulation (e.g. "data")


    """
    def __init__(self, filename, delimiter="\t", output_dir="."):

        self.filename = filename
        self.matrix = []
        self.header = []

        row_num = 0
        with open(self.filename, newline='') as file:
            reader = csv.reader(file, delimiter='\t')
            self.header = reader.__next__()

            for row in reader:
                row_num += 1

                # data check: all rows should have same length as the header row
                if not len(row) == len(self.header):
                    print("err:", len(row), str(row_num), row)

                # Unicode normalization (form D) and strip whitespace
                for token in row:
                    token = unicodedata.normalize("NFD", token)
                    token = token.strip()

                self.matrix.append(row)


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

    def pprint(self, matrix):
        """
        Print the matrix in tab delimited format
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


if __name__=="__main__":
    s = Spreadsheet("/Users/stiv/Dropbox/Projects/dogon/Moran_Dogon.comp.vocab.UNICODE.csv")
    # s = Spreadsheet("/Users/stiv/Projects/lingpy/test/spreadsheet.tsv")

    # s.print_matrix_stats()
    s.print_qlc_format()

    # frm = s.get_matrix_full_rows()
    # s.pprint(frm)
    # print(len(frm))
    # s.print_qlc_format(frm)

