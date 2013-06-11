# -*- coding: utf-8 -*-

"""
Work in progress of a module to handle some of the tokenization of data.

Will be superseded.

@date: 2012-06-01
@author: Steven Moran

"""

import os
import sys
import unicodedata
import collections

from qlc.corpusreader import CorpusReaderWordlist
from qlc.orthography import OrthographyParser
from qlc.orthography import GraphemeParser
from configparser import SafeConfigParser
from configparser import ConfigParser

# class tokenzier:
# read in config file
# do data manipulation / storage
# returns various data types: list of concepts tokenized...
# use to produce various things, like initial orthography profile


# NOTE: if you have more than one class at the same level in Python; both are fired if you create an object in main.

class SpreadsheetParser:
    pass

    def __init__(self): 
        input_file = open("data/dogon/heath2012_spreadsheet_format.csv", "r")
        # default read header - 1st row is 

    # size of spreadsheet (x*y)

    # number of filled cells total

    # number of cells per resource (col1 should equal number of rows, i.e. all should be filled)

    """
    output_file = open("data/dogon/qlc_output_heath2012.tsv", "w")
    header = input_file.readline()
    header = header.strip()
    languages = header.split("\t")

    # print("COUNTERPART"+"\t"+"CONCEPT"+"\t"+"LANGUAGE_BOOKNAME")
    # actually: concept, counterpart, language
    
    # this is essentially Mattis's format
    for line in input_file:
        line = line.strip()
        line = unicodedata.normalize("NFD", line) # skip the END
        tokens = line.split("\t")
        concept = tokens[0].strip()
        for i in range(1, len(tokens)-1):
            counterpart = tokens[i].strip()
            if counterpart == "":
                counterpart = "NONE"
            result = counterpart+"\t"+concept+"\t"+languages[i].strip()
            # result = concept+"\t"+tokens[i].strip()+"\t"+languages[i].strip()
            # print(result)
            output_file.write(result+"\n")
    output_file.close()

    input_file.close()
    """

    

class Tokenizer:

    """ takes as input a file with the QLC format:
    counterpart \t concept \t language

    and does things like 

    - tokenizes the file into LINGPY format
    - tokenizes the data into ortographically parsed QLC format
    - locates unicorns

    """


    def __init__(self):
        # deal with configuration file
        # configparser.read(default.cfg)
        cfg = SafeConfigParser()
        cfg.read("default.cfg")

        data = cfg.get("Paths", "data")
        orthography_profile = cfg.get("Paths", "orthography_profile")

        # set variables, e.g. source, orthography parser, etc.
        self.data = open(data, "r")

        self.o = OrthographyParser(orthography_profile)        
        # self.o = GraphemeParser()        

        self._languages = collections.defaultdict(int) # given unique ID to each unique language name
        self._concepts = collections.defaultdict(int) # ...
        self._counterparts = collections.defaultdict(int) # ..
        self._wordlist_iterator = self._process_input(self.data)

        # print(type(self.iterator))
        # print(len(self.counterparts))
        # words = self.get_qlc_tokenized_words()

        """
        count = 0
        for line in words:
            if line != "":
                print(line)
                count += 1
        print(count)
        """

        """
        self.cr = CorpusReaderWordlist("data/csv")
        self.wordlist_iterator = ( (wordlistdata_id, concept, counterpart)
            for wordlistdata_id in self.cr.wordlistdata_ids_for_bibtex_key(source)
            for concept, counterpart in self.cr.concepts_with_counterparts_for_wordlistdata_id(wordlistdata_id)
        )
        """

    def _process_input(self, file):
        languages_id = 1
        concepts_id = 1
        counterparts_id = 1
        header = file.readline()

        lines = []
        for line in file:
            line = line.strip()
            line = line.replace("  ", " ")
            counterpart, concept, language = line.split("\t")
            result = (counterpart, concept, language)
            lines.append(result)

            if language not in self._languages:
                self._languages[language] = languages_id
                languages_id += 1
            if concept not in self._concepts:
                self._concepts[concept] = concepts_id
                concepts_id += 1
            if counterpart not in self._counterparts:
                counterparts_id += 1
                self._counterparts[counterpart] = counterparts_id

        return ((concept, counterpart, language) for concept, counterpart, language in lines)


    def get_qlc_tokenized_words(self):
        unparasables = open("unparsables.txt", "w")
        tokenized_words = []
        for counterpart, concept, language in self._wordlist_iterator:
            counterpart = unicodedata.normalize("NFD", counterpart)
            grapheme_parsed_counterpart_tuple = self.o.parse_string_to_graphemes_string(counterpart)
            if grapheme_parsed_counterpart_tuple[0] == False:
                unparsables.write(grapheme_parsed_counterpart_tuple[1])
                continue
        
            grapheme_parse = grapheme_parsed_counterpart_tuple[1]
            tokenized_words.append(grapheme_parse)
        return tokenized_words

    def get_ipa_tokenized_words(self):
        tokenized_words = []
        words = get_list_qlc_tokenized_words()
        for word in words:
            grapheme_parsed_counterpart_tuple = self.o.parse_string_to_graphemes_string(counterpart)
            
    def lingpy_output(self):
        row_id = 1
        # given some data set from the corpusreader or somewhere else, output a lingpy format
        print("ID"+"\t"+"Taxa"+"\t"+"TaxonID"+"\t"+"Gloss"+"\t"+"GlossID"+"\t"+"IPA"+"\t"+"Orthography")
        # print("# LANGUAGE"+"\t"+"CONCEPT"+"\t"+"COUNTERPART"+"\t"+"ORTHO_PARSE")

        for counterpart, concept, language in self._wordlist_iterator:
            # counterpart, concept, language in self._wordlist_iterator:
            # skip for Mattis
            if counterpart == "?" or counterpart == "NONE":
                continue

            grapheme_parsed_counterpart_tuple = self.o.parse_string_to_graphemes_string(counterpart)
            if grapheme_parsed_counterpart_tuple[0] == False:
                continue

            ortho_parse = grapheme_parsed_counterpart_tuple[1]

            print(str(row_id)+"\t"+language+"\t"+str(self._languages[language])+"\t"+concept+"\t"+str(self._concepts[concept])+"\t"+counterpart+"\t"+grapheme_parsed_counterpart_tuple[1])
            # print(language+"\t"+concept+"\t"+counterpart+"\t"+grapheme_parsed_counterpart_tuple[1])

            row_id += 1

    def matrix_output(self):
        # produce Jelena style output format with matrix
        pass

    def qlc_output_format(self):
        # produce counterpart \t concept \t language QLC output format
        print("COUNTERPART"+"\t"+"ORTHO_PARSE"+"\t"+"CONCEPT"+"\t"+"LANGUAGE")
        for counterpart, concept, language in self._wordlist_iterator:
            if counterpart == "?" or counterpart == "NONE":
                print(counterpart+"\t"+counterpart+"\t"+concept+"\t"+language)                
                continue
            grapheme_parsed_counterpart_tuple = self.o.parse_string_to_graphemes_string(counterpart)
            
            # skip shit that doesn't parse
            if grapheme_parsed_counterpart_tuple[0] == False:
                continue

            ortho_parse = grapheme_parsed_counterpart_tuple[1]
            print(counterpart+"\t"+ortho_parse+"\t"+concept+"\t"+language)




if __name__=="__main__":
    from qlc.tokenizer import Tokenizer
    from qlc import ngram
    t = Tokenizer()
    t.lingpy_output()
#    words = t.get_qlc_tokenized_words()
#    ngram.unigram_model(words)
