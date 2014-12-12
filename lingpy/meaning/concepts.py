# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

# author   : Peter Bouda
# email    : pbouda@cidles.eu
# created  : 2013-06-06 11:57
"""
This module provides basic classes and functions for the handling of concepts.

TODO:
* How to handle concept comparison in ConceptGraph?

"""

__author__="Peter Bouda"
__date__="2013-07-29"

import os
import sys
from datetime import date,datetime
import numpy as np
import pickle
import codecs
import re
from operator import itemgetter
import abc
import csv
import glob
import zipfile

# basic lingpy imports
from ..settings import rcParams
from ..basic.dictionary import Dictionary
from .. import log
from .. import util

try:
    from nltk.stem.snowball import SpanishStemmer
except ImportError:
    log.missing_module('nltk')


class ConceptGraph():
    """
    Basic class to represent concepts and lexical items in a network.

    Parameters
    ----------
    concepts : list
        A list of concepts to initialize the network. The type of the concecpts
        is str, in the standard case, but concepts might be encoded by
        any type that can be used as a dictionary key. In any case the user
        has to make sure that the ConceptComparer class uses the same type for
        semantic comparison.
    pivot_lang_iso : str
        The ISO 639-2 code of the pivot language. The code is used when the user
        adds a new dictionary, to extract the correct entries for comparison
        from the dictionary.
    concept_comparer: ConceptComparerBase
        A sub-class of ConceptComparerBase, that is used to compare lexical
        items against the concepts.

    Notes
    -----
    The ConceptGraph class is used to represent a network of concepts and
    lexical items. At the moment it is used to connect a number of
    :py:class:`~lingpy.lingpy.basic.dictionary.Dictionary`
    via their translation. For example, if two or more dictionaries share one
    language on the head or translation side (e.g. all have translations in
    Spanish), the entries can be linked via this pivot language.

    Initially the ConceptGraph consists only of nodes for the concepts (for
    example the concept of a Swadesh list). Then, lexical items can be linked
    to those concepts based on semantic comparison. To carry out the semantic
    comparison you will use a sub-class of
    :py:class:`~lingpy.meaning.concepts.ConceptComparerBase`. Lingpy currently
    provides one implementation for semantic comparison based on the comparison
    of Spanish word stems agains the stems of the Spanish Swadesh List. The
    class :py:class:`~lingpy.meaning.concepts.ConceptComparerSpanishStem`
    provides this implentation.

    When all lexical items are linked against the concept the result can be
    exported as a lingpy :py:class:`~lingpy.lingpy.basic.wordlist.Wordlist`.

    Examples
    --------
    See script in `scripts/dictionary/extract_component_as_wordlist.py`, which
    extracts a wordlist from several qlc dictionaries. This is a nice
    demonstration of the whole workflow.

    See also
    --------
    ConceptComparerBase
    ConceptComparerSpanishStem
    lingpy.basic.dictionary.Dictionary
    lingpy.basic.wordlist.Wordlist

    """

    def __init__(self, concepts, pivot_lang_iso, concept_comparer):
        self.pivot_lang_iso = pivot_lang_iso
        self.concept_comparer = concept_comparer
        self.graph = {}
        self.use_tokens = False
        self.doculects = set()
        for concept in concepts:
            self.graph[concept] = set()

    def add_dictionary(self, dictionary):
        """
        Adds the lexical items of a dictionary to the concept graph.


        Parameters
        ----------
        dictionary : lingpy.basic.dictionary.Dictionary
            The dictionary to add to the concept graph.

        Returns
        -------
        Nothing.

        Notes
        -----
        Add the entries of a
        :py:class:`~lingpy.lingpy.basic.dictionary.Dictionary` to the concept
        graph. The methods checks which part of the dictionary matches
        the pivot language. If none of the languages matches, then no
        items are added to the concept graph.

        If the dictionary was already tokenized with an orthography profile,
        then those tokens will be added to the concept graph and linked to the
        concepts. Otherwise plain strings are used, either the heads or
        translations, depending on the pivot language.

        To compare the lexical items to the concepts the concept comparer that
        was used to initialize the class is called (a sub-class of
        :py:class:`~lingpy.meaning.concepts.ConceptComparerBase`).

        See also
        --------
        ConceptComparerBase
        lingpy.basic.dictionary.Dictionary

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
                if self.concept_comparer.compare_to_concept(pivot, concept):
                    self.graph[concept].add((entry[0], trans, tokens, doculect))

        for doculect, iso in dictionary.doculect2iso.items():
            self.doculects.add((doculect, iso))

    def output_wordlist(self, stream):
        """
        Export a wordlist from the concept graph.

        Parameters
        ----------
        filename : IO stream or str
            A path to the filename to write the wordlist. The wordlist is
            written as an qlc wordlist.

        Returns
        -------
        Nothing.

        Notes
        -----
        The wordlist can be read in again with the class
        :py:class:`~lingpy.lingpy.basic.wordlist.Wordlist`. If the dictionaries
        that were added already contained tokenized string (for example with
        support of an orthography profile), then the resulting wordlist will
        also contain a column "tokes" with the tokenized data.

        See also
        --------
        lingpy.basic.wordlist.Wordlist

        """
        if not hasattr(stream, 'read'):
            stream = codecs.open(stream, "w", "utf-8")

        # write header
        stream.write("@date: {0}\n".format(str(date.today())))

        stream.write(
            "@source_title: Automatically created wordlist, by lingpy.\n")

        for doculect, iso in self.doculects:
            stream.write("@doculect: {0}, {1}\n".format(doculect, iso))

        if self.use_tokens:
            stream.write(
                "QLCID\tCONCEPT\tCOUNTERPART\tCOUNTERPART_DOCULECT\tTOKENS\n")
            for concept in self.graph:
                    for qlcid, counterpart, tokens, counterpart_doculect \
                            in self.graph[concept]:
                        stream.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                            qlcid, concept.upper(), counterpart,
                            counterpart_doculect, tokens))
        else:
            stream.write(
                "QLCID\tCONCEPT\tCOUNTERPART\tCOUNTERPART_DOCULECT\n")

            for concept in self.graph:
                for qlcid, counterpart, _, counterpart_doculect \
                        in self.graph[concept]:
                    stream.write("{0}\t{1}\t{2}\t{3}\n".format(
                        qlcid, concept.upper(), counterpart,
                        counterpart_doculect))

        stream.close()


class ConceptComparerBase():
    """ 
    Base class for semantic comparison. 

    Notes
    ------
    Create a sub-class of this if you want
    to implement your own methods to compare lexical items against concepts.

    The base class forces the implementation of one method
    `compare_to_concept()`. This method is responsible to compare a given
    lexical entry (for example as string) against a given concept (that might
    be a string or any other type).

    Lingpy currently implements one concept comparer that compares a given
    lexical items against a given concept based on a spanish stemmer as
    class :py:class:`~lingpy.meaning.concepts.ConceptComparerSpanishStem`.

    The class is used for example by the class 
    :py:class:`~lingpy.meaning.concepts.ConceptGraph` to link concept and
    lexical items.

    See also
    --------
    ConceptComparerSpanishStem
    ConceptGraph
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def compare_to_concept(self, element, concept):
        """Compares a given element to a concept. Implement this for your
        type of concept.

        Parameters
        ----------
        element : str
            The string (for example a lexical item: head or translation) to
            compare to the concept.
        concept : str or object
            The conpect to compare to.

        Return
        ------
        match : bool
            True if element matches the given concept, False otherwise.

        """
        raise NotImplementedError("Method must be implemented")

class ConceptComparerStringMatch(ConceptComparerBase):
    """
    Implementation of a concept comparer based on a stemmer for spanish.

    Parameters
    ----------
    None.

    Notes
    -----
    This is a sub-class of
    :py:class:`~lingpy.meaning.concepts.ConceptComparerBase`. It uses a simple
    match of the stem of a given (spanish) string against a given context
    (that is supposed to be a stemmed spanish word stem).

    See also
    --------
    ConceptComparerBase
    ConceptGraph
    """

    def __init__(self, complete=True):
        self.re_brackets = re.compile(" ?\([^)]\)")
        self.complete = complete

    def compare_to_concept(self, element, concept):
        """Compares a given element to a concept.

        Parameters
        ----------
        element : str
            The string (for example a lexical item: head or translation) to
            compare to the concept.
        concept : str or object
            The conpect to compare to.

        Return
        ------
        match : bool
            True if element matches the given concept, False otherwise.

        Notes
        -----
        The `element` is supposed to be a spanish word, the concept a string,
        for example from a Swadesh list.

        See also
        --------
        spanish_swadesh_list

        """
        element = self.re_brackets.sub("", element)
        element = element.strip()
        if self.complete:
            if concept == element:
                return True
        else:
            print(element)
            if concept in element:
                return True
        return False

class ConceptComparerSpanishStem(ConceptComparerBase):
    """
    Implementation of a concept comparer based on a stemmer for spanish.

    Parameters
    ----------
    None.

    Notes
    -----
    This is a sub-class of
    :py:class:`~lingpy.meaning.concepts.ConceptComparerBase`. It uses a simple
    match of the stem of a given (spanish) string against a given context
    (that is supposed to be a stemmed spanish word stem).

    See also
    --------
    ConceptComparerBase
    ConceptGraph
    """

    def __init__(self):
        self.stemmer = SpanishStemmer(True)
        self.re_brackets = re.compile(" ?\([^)]\)")

    def compare_to_concept(self, element, concept):
        """Compares a given element to a concept.

        Parameters
        ----------
        element : str
            The string (for example a lexical item: head or translation) to
            compare to the concept.
        concept : str or object
            The conpect to compare to.

        Return
        ------
        match : bool
            True if element matches the given concept, False otherwise.

        Notes
        -----
        The `element` is supposed to be a spanish word, the concept a stemmed
        entry of the spanish Swadesh List.

        See also
        --------
        spanish_swadesh_list

        """
        element = self.re_brackets.sub("", element)
        element = element.strip()
        if not " " in element:
            stem = self.stemmer.stem(element)
            if stem == concept:
                return True
        return False

def spanish_swadesh_list(stemmed=True):
    """
    Helper function that returns a list of strings with the stems of the
    spanish Swadesh entries.

    """
    try:
        stemmer = SpanishStemmer(True)
    except:
        log.warn("Spanish stemmer could not be loaded!")
        return

    swadesh_entries = []
    for line in util.read_text_file(
            util.data_path('swadesh', 'swadesh_spa.txt'), lines=True):
        line = line.strip()
        for e in line.split(","):
            e = e.strip()
            if stemmed:
                stem = stemmer.stem(e)
                swadesh_entries.append(stem)
            else:
                swadesh_entries.append(e)
    return swadesh_entries


def extract_component_as_wordlist(component, use_profiles, concepts, iso_pivot,
    output_file=None, complete=True):
    """
    This function extracts a wordlist from the quanthistling dictionaries.
    The package of dictionaries should be extracted to the current working
    directory first. If no "sources.csv" file is found the function will try
    to download the filtered data package, i.e. only the entries that contain
    a word from the Spanish Swadesh list will be used.

    The resulting wordlist will be written to a file.

    Parameters
    ----------
    component : str or list
        The component for which to extract the dictionaries. "Witotoan", for
        example. This can also be a list of all bibtex_keys of the books today
        be processed.
    use_profiles : bool
        Whether to apply orthography profiles or not. If true, than any
        dictionary for which there is no profile in lingpy will be skipped.
    concepts : list
        A list of concepts to initialize the network. The type of the concecpts
        is str.
    iso_pivot : str
        The ISO 639-2 code of the pivot language. The code is used when the user
        adds a new dictionary, to extract the correct entries for comparison
        from the dictionary.
    output_file: str
        The filename of the wordlist to output. If none is given, then the
        component will be used for the filename (e.g. "Witototan.csv").
    complete : bool
        Match concepts by exact match (e.g. "fuego" == "fuego") if True. If
        False then the concepts might also be part of the entry (e.g. "fuego"
            == "asar calentando algo encima del fuego a medio cocer.")

    Return
    ------
    output_file : str
        The filename of the output file.

    """
    if type(component) is list:
        component_sources = component
    else:
        if not os.path.exists("sources.csv"):
            try:
                import requests
            except:
                print("Module 'requests' not found. I cannot download the data "
                      "automatically for you.\n\nPlease download manually at:\n"
                      "http://www.quanthistling.info/data/downloads/csv/data.zip")
                return False

            r = requests.get(
                "http://www.quanthistling.info/data/downloads/csv/data.zip")
            with open("data.zip", "wb") as f:
                f.write(r.content)

            z = zipfile.ZipFile("data.zip")
            z.extractall()

        sources = csv.reader(codecs.open("sources.csv", "r", "utf-8"),
            delimiter="\t")
        component_sources = list()
        for source in sources:
            # add for "ready"-only dicts: and source[3] == "True"
            if source[5] == component and source[1] == "dictionary":
                component_sources.append(source[0])

    cm = ConceptComparerStringMatch(complete)
    cg = ConceptGraph(concepts, iso_pivot, cm)

    for f in glob.glob("*.csv"):
        if ("-" in f and f[:f.index("-")] in component_sources) or \
                ("." in f and f[:f.index(".")] in component_sources):
            print("Adding {0}...".format(f))
            di = Dictionary(f)
            if use_profiles:
                ortho_path = util.data_path(
                    'orthography_profiles', "{0}.prf".format(f[:f.index("-")]))
                if os.path.exists(ortho_path):
                    if iso_pivot in di.head_iso:
                        di.tokenize(ortho_path, source="translation")
                    else:
                        di.tokenize(ortho_path)
                    cg.add_dictionary(di)
                else:
                    print(
                        "  Orthography profile not found, skipping dictionary.")
            else:
                cg.add_dictionary(di)

    if not output_file:
        output_file = "{0}.csv".format(component)
    cg.output_wordlist(output_file)
    return output_file
