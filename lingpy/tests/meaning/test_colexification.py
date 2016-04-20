"""
Tests for colexification module.
"""
from six import text_type
from lingpy.tests.util import WithTempDir 
from lingpy.tests.util import test_data
from lingpy.meaning.colexification import *
import lingpy.meaning.colexification as colx
from lingpy.basic.wordlist import Wordlist


class TestColexifications(WithTempDir):
    
    def setUp(self):
        WithTempDir.setUp(self)
        self.wordlist = Wordlist(test_data('colexification.tsv'))
        self.cols = colx._get_colexifications(self.wordlist)
    
    def test_colexification_network(self):
        graph = colexification_network(Wordlist(test_data('colexification.tsv')))
        assert "hand" in graph and "arm" in graph
        
        graph = colexification_network(Wordlist(test_data('colexification.tsv')),
                bipartite=True)
        assert 'arm' in graph.edge['l4.4'] and 'hand' in graph.edge['l4.4']
        
        graph = colexification_network(Wordlist(test_data('colexification.tsv')),
                output="gml", filename=text_type(self.tmp_path("test")))
        
    
    def test__get_colexifications(self):
        
        assert len(self.cols[0]) == 5
    
    def test__get_colexifications_by_taxa(self):
        
        colt = colx._get_colexifications_by_taxa(self.cols)
        assert ('arm', 'hand') in colt['l1']
    
    def test__get_statistics(self):
        
        stats = colx._get_statistics(self.wordlist)
        assert stats['arm']['wordOcc'] == 6
    
    def test__make_matrix(self):
        
        colt = colx._get_colexifications_by_taxa(self.cols)
        taxa = ['l1', 'l2', 'l3', 'l4', 'l5', 'l6']
        
        matrix = colx._make_matrix(taxa, colt)
        assert matrix[0][0] == 0
    
    def test__make_graph(self):
        graph = colx._make_graph(self.cols, bipartite=True)
        assert 'arm' in graph.edge['l4.4'] and 'hand' in graph.edge['l4.4']

    def test_compare_colexifications(self):

        matrix = colx.compare_colexifications(self.wordlist)
        assert matrix[0][0] == 0

    def test_evaluate_colexifications(self):
        graph = colx._make_graph(self.cols)
        nodes, edges = colx.evaluate_colexifications(graph, weight='wordWeight')
        colx.evaluate_colexifications(graph, weight='wordWeight',
                outfile=text_type(self.tmp_path('test')))


