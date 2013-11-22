# author   : Peter Bouda
# email    : pbouda@cidles.eu
# created  : 2013-08-26 09:48
"""
This contains the test classes and functions for concepts.py.

"""

__author__="Peter Bouda"
__date__="2013-08-26"

import os
import tempfile
import filecmp

import lingpy
import lingpy.meaning.concepts

class TestConceptComparerSpanishStem:

    def setup(self):
    	self.cm = lingpy.meaning.concepts.ConceptComparerSpanishStem()

    def test_compare_to_concept(self):
        res = self.cm.compare_to_concept("comer", "com")
        assert res == True

        res = self.cm.compare_to_concept("viver", "cas")
        assert res == False


class TestConceptGraph:

    def setup(self):
        concepts = lingpy.meaning.concepts.spanish_swadesh_list()
        cm = lingpy.meaning.concepts.ConceptComparerSpanishStem()
        
        self.cg = lingpy.meaning.concepts.ConceptGraph(concepts, "spa", cm)

    def test_add_dictionary(self):
        dictionary = lingpy.Dictionary(os.path.join(os.path.dirname( __file__ ),
            '..', 'test_data', 'burtch1983-19-262.csv'))
        self.cg.add_dictionary(dictionary)

        assert len(self.cg.graph) == 210

        assert self.cg.doculects == {('Castellano', 'spa'),
            ('Huitoto Murui', 'huu')}

    def test_output_wordlist(self):
        dictionary = lingpy.Dictionary(os.path.join(os.path.dirname( __file__ ),
            '..', 'test_data', 'burtch1983-19-262.csv'))
        self.cg.add_dictionary(dictionary)

        dictionary = lingpy.Dictionary(os.path.join(os.path.dirname( __file__ ),
            '..', 'test_data', 'aleman2000-9-69.csv'))
        self.cg.add_dictionary(dictionary)

        assert len(self.cg.graph) == 210

        assert self.cg.doculects == {('Desano', 'des'), ('Español', 'spa'),
            ('Castellano', 'spa'), ('Huitoto Murui', 'huu')}

        _, tmp_filename = tempfile.mkstemp()
        wordlist = os.path.join(os.path.dirname( __file__ ),
            '..', 'test_data', 'concept_graph_wordlist_out.csv')
        self.cg.output_wordlist(tmp_filename)

        filecmp.cmp(tmp_filename, wordlist)

