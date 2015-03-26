# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
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

from six import text_type as str

# add a warning routine to prevent this from firing
    
class TestConceptComparerSpanishStem:

    def setup(self):
        try:
    	    self.cm = lingpy.meaning.concepts.ConceptComparerSpanishStem()
        except:
            self.cm = None

    def test_compare_to_concept(self):

        if self.cm:
            res = self.cm.compare_to_concept("comer", "com")
            assert res == True

            res = self.cm.compare_to_concept("viver", "cas")
            assert res == False

class TestConceptComparerStringMatch:

    def setup(self):
        try:
            self.cm = lingpy.meaning.concepts.ConceptComparerStringMatch()
        except:
            self.cm = None
    def test_compare_to_concept(self):
        if self.cm:
            res = self.cm.compare_to_concept("comer", "comer")
            assert res == True
            res = self.cm.compare_to_concept("Tu comer", "comer")
            assert res == False

            self.cm.complete = False
            res = self.cm.compare_to_concept("comer", "comer")
            assert res == True
            res = self.cm.compare_to_concept("Tu comer", "comer")
            assert res == True
            self.cm.complete = True

            res = self.cm.compare_to_concept("viver", "cas")
            assert res == False

class TestConceptGraph:

    def setup(self):
        try:
            concepts = lingpy.meaning.concepts.spanish_swadesh_list()
            cm = lingpy.meaning.concepts.ConceptComparerSpanishStem()
        
            self.cg = lingpy.meaning.concepts.ConceptGraph(concepts, "spa", cm)
        except:
            self.cg = None

    def test_add_dictionary(self):
        if self.cg:
            dictionary = lingpy.Dictionary(os.path.join(os.path.dirname( __file__ ),
                '..', 'test_data', 'burtch1983-19-262.csv'))
            self.cg.add_dictionary(dictionary)

            assert len(self.cg.graph) == 210

            assert self.cg.doculects == {('Castellano', 'spa'),
                ('Huitoto Murui', 'huu')}

    def test_output_wordlist(self):
        if self.cg:
            dictionary = lingpy.Dictionary(os.path.join(os.path.dirname( __file__ ),
                '..', 'test_data', 'burtch1983-19-262.csv'))
            self.cg.add_dictionary(dictionary)

            dictionary = lingpy.Dictionary(os.path.join(os.path.dirname( __file__ ),
                '..', 'test_data', 'aleman2000-9-69.csv'))
            self.cg.add_dictionary(dictionary)

            assert len(self.cg.graph) == 210
            
            assert ('Desano','des') in self.cg.doculects
            assert ('Castellano','spa') in self.cg.doculects

            _, tmp_filename = tempfile.mkstemp()
            wordlist = os.path.join(os.path.dirname( __file__ ),
                '..', 'test_data', 'concept_graph_wordlist_out.csv')
            self.cg.output_wordlist(tmp_filename)

            filecmp.cmp(tmp_filename, wordlist)

