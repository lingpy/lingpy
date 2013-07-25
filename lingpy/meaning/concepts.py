# author   : Peter Bouda
# email    : pbouda@cidles.eu
# created  : 2013-06-06 11:57
"""
This module provides basic classes and functions for the handling of concepts.

TODO:
* How to handle concept comparison in ConceptGraph?

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
from ..settings import rcParams

try:
    from nltk.stem.snowball import SpanishStemmer
except ImportError:
   print(rcParams['W_missing_module'].format("nltk"))

class ConceptGraph():
    """

    """

    def __init__(self, concepts, pivot_lang_iso, concept_matcher):
        """

        """
        #self.concepts = concepts
        self.pivot_lang_iso = pivot_lang_iso
        self.concept_matcher = concept_matcher
        self.graph = {}
        self.use_tokens = False
        self.doculects = set()
        for concept in concepts:
            self.graph[concept] = set()

    def add_dictionary(self, dictionary):
        """

        """
        columns = [ "qlcid", "head", "translation", "head_doculect",
                    "translation_doculect" ]
        if "tokens" in dictionary.header:
            columns.append("tokens")
            self.use_tokens = True
        for entry in dictionary.get_tuples(columns):
            pivot = ""; trans = ""; doculect = ""
            if self.pivot_lang_iso in dictionary.head_iso:
                pivot = entry[1]
                trans = entry[2]
                doculect = entry[4]
            elif self.pivot_lang_iso in dictionary.translation_iso:
                pivot = entry[2]
                trans = entry[1]
                doculect = entry[3]
            else:
                continue

            tokens = ""
            if self.use_tokens:
                tokens = ' '.join(entry[5])

            for concept in self.graph:
                if self.concept_matcher.compare_to_concept(pivot, concept):
                    self.graph[concept].add((entry[0], trans, tokens, doculect))

        for doculect, iso in dictionary.doculect2iso.items():
            self.doculects.add((doculect, iso))

    def output_wordlist(self, filename):
        """

        """
        wordlist = codecs.open(filename, "w", "utf-8")

        # write header
        wordlist.write("@date: {0}\n".format(str(date.today())))

        wordlist.write(
            "@source_title: Automatically created wordlist, by lingpy.\n")

        for doculect, iso in self.doculects:
            wordlist.write("@doculect: {0}, {1}\n".format(doculect, iso))

        if self.use_tokens:
            wordlist.write(
                "QLCID\tCONCEPT\tCOUNTERPART\tCOUNTERPART_DOCULECT\tTOKENS\n")
            for concept in self.graph:
                    for qlcid, counterpart, tokens, counterpart_doculect \
                            in self.graph[concept]:
                        wordlist.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                            qlcid, concept.upper(), counterpart,
                            counterpart_doculect, tokens))
        else:
            wordlist.write(
                "QLCID\tCONCEPT\tCOUNTERPART\tCOUNTERPART_DOCULECT\n")

            for concept in self.graph:
                for qlcid, counterpart, _, counterpart_doculect \
                        in self.graph[concept]:
                    wordlist.write("{0}\t{1}\t{2}\t{3}\n".format(
                        qlcid, concept.upper(), counterpart,
                        counterpart_doculect))

        wordlist.close()


class ConceptComparerBase():
    """
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def compare_to_concept(self, element, concept):
        """Compares a given element to a concept. Implement this for your
        type of concept.

        Parameters
        ----------
        element : str
            The string to compare to the concept.
        concept : str or object
            The conpect to compare to.

        Return
        ------
        match : bool
            True if element matches the given concept, False otherwise.

        """
        raise NotImplementedError("Method must be implemented")

class ConceptComparerSpanishStem(ConceptComparerBase):

    def __init__(self):
        self.stemmer = SpanishStemmer(True)
        self.re_brackets = re.compile(" ?\([^)]\)")

    def compare_to_concept(self, element, concept):
        element = self.re_brackets.sub("", element)
        element = element.strip()
        if not " " in element:
            stem = self.stemmer.stem(element)
            if stem == concept:
                return True
        return False

def spanish_swadesh_list():
    stemmer = SpanishStemmer(True)
    # load swadesh list
    swadesh_file = os.path.split(
                    os.path.dirname(
                        os.path.abspath(
                            __file__
                            )
                        )
                    )[0] + '/data/swadesh/swadesh_spa.txt'

    swadesh = codecs.open(swadesh_file, "r", "utf-8")

    swadesh_entries = []
    for line in swadesh:
        line = line.strip()
        for e in line.split(","):
            e = e.strip()
            stem = stemmer.stem(e)
            swadesh_entries.append(stem)
    return swadesh_entries
