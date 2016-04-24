from lingpy.tests.util import test_data
from lingpy.read.csv import *
from lingpy.meaning.glosses import parse_gloss, compare_conceptlists,\
        compare_concepts
from lingpy.tests.util import WithTempDir
from six import text_type

class Tests(WithTempDir):

    def test_parse_gloss(self):
        data = csv2list(test_data('Ringe-2002-421.tsv'))[1:]
        target_list = [x[0] for x in csv2list(test_data('target_list_ringe.tsv'))]
        for line, target in zip(data, target_list):
            datum = line[1]
            glosses = parse_gloss(datum)
            for a, b, c, d, e, f, g, h, i in glosses:
                # print(datum,'=>',','.join([x for x in [a,b,c,d,e,f,g,''.join(h),i]]),'\t',a)
                assert a == target
    
    
    def test_conceptlists(self):
        comparison1 = compare_conceptlists(test_data('listA.tsv'), test_data('listB.tsv'))
        comparison2 = compare_conceptlists(test_data('listA.tsv'),
                test_data('listB.tsv'), output="tsv",
                filename=text_type(self.tmp_path('out')))
        
        assert comparison2 is None
        assert isinstance(comparison1, list)
    
    def test_compare_concepts(self):

        c1 = 'hand (noun)'
        c2 = 'hand (or arm)'
        c3 = 'hand (noun, also arm)'
        c4 = 'hand / arm'
        c5 = 'hand arm'
        c6 = 'hand and arm (noun)'
        c7 = 'the hand or arm'
        
        sims1 = compare_concepts(c1, c2)
        sims2 = compare_concepts(c1, c1)
        sims3 = compare_concepts(c1, c3)
        sims4 = compare_concepts(c4, c5)
        sims5 = compare_concepts(c6, c7)
        sims6 = compare_concepts(c5, c4)

        
        assert sims1[0][-1] == 3
        assert sims2[0][-1] == 1
        assert sims3[0][-1] == 2
        assert sims4[0][-1] == 5
        assert sims5[0][-1] == 4
        assert sims6[1][-1] == 6

        

