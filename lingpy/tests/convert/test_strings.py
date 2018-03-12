"""
Test conversions involving strings.
"""
from __future__ import division, unicode_literals
from unittest import TestCase
from six import text_type
from lingpy.tests.util_testing import WithTempDir

import sys
import lingpy
from lingpy.basic.wordlist import Wordlist
from lingpy.convert.strings import scorer2str, msa2str, matrix2dst, pap2nex, pap2csv, write_nexus, _is_constant
from lingpy.util import read_text_file
from lingpy.tests.util import test_data


class Tests(TestCase):
    def test_scorer2str(self):
        """
        Test conversion of scorers to strings.
        """
        self.assertEqual(scorer2str(lingpy.rc('dolgo').scorer),
                         read_text_file(test_data('dolgo.scorer')))

    def test_msa2str(self):
        aranger = '{body}{meta}'

        # read msa traditionally into an object
        msaA = lingpy.MSA(test_data('harry.msa'))

        # read msa from dictionary
        msaB = lingpy.read.qlc.read_msa(test_data('harry.msa'))

        # read msa with IDs
        msaC = lingpy.read.qlc.read_msa(
            test_data('harry_with_ids.msa'), ids=True, header=False)

        # we adjust the dataset and the seq_id since otherwise we won't have
        # similar output
        msaC['seq_id'] = 'test'
        msaC['dataset'] = 'file'

        # when converting these different objects to string with the same body and
        # the like, they should be identical, so we check this here
        strA = msa2str(msaA, _arange=aranger)
        strB = msa2str(msaB, _arange=aranger)
        strC = msa2str(msaC, _arange=aranger, wordlist=False)

        assert strA == strB == strC

        # we next test for converting with the merging attribute
        strD = msa2str(msaC, _arange=aranger, wordlist=True, merge=True)
        strE = msa2str(msaC, _arange=aranger, wordlist=True, merge=False)

        # remove tabstops for checking similar strings
        strDst = strD.replace('\t', '')
        strEst = strE.replace('\t', '')

        # get index until 'COLUMN'
        idx = strDst.index('COLUMNID')
        assert strD != strE and strDst[:idx] == strEst[:idx]

        # add a consensus string to all msa objects
        consensusA = lingpy.align.sca.get_consensus(lingpy.align.sca.MSA(msaB), gaps=True)
        consensusB = lingpy.align.sca.get_consensus(lingpy.align.sca.MSA(msaC), gaps=True)

        msaB['consensus'] = consensusA
        msaC['consensus'] = consensusB

        assert msa2str(msaB) == msa2str(msaC, wordlist=False)

    def test_matrix2dst(self):
        matrix = lingpy.algorithm.squareform([0.5, 0.75, 0.8])

        # we choose same format for taxa as default
        taxa = ['t_1', 't_2', 't_3']

        phylA = matrix2dst(matrix, taxa=taxa)
        phylB = matrix2dst(matrix)

        assert phylA == phylB

        phylC = matrix2dst(matrix, taxa=taxa, stamp='# Written with joy.')
        phylD = matrix2dst(matrix, stamp='# Written with joy.')

        assert phylC == phylD

        phylE = matrix2dst(matrix, taxa=taxa, taxlen=20)
        phylF = matrix2dst(matrix, taxlen=30)

        assert 18 * ' ' in phylE and 28 * ' ' in phylF

        # check for tab-stop output when taxlen is set to 0
        self.assertEqual(matrix2dst(matrix, taxlen=0).count('\t'), 9)

    def test_pap2nex(self):
        nex = """#NEXUS

BEGIN DATA;
DIMENSIONS ntax=2 NCHAR=4;
FORMAT DATATYPE=STANDARD GAP=- MISSING=0 interleave=yes;
MATRIX

a 1111
b 0101

;

END;"""

        taxa = ['a', 'b']
        papsA = [
            [1, 0],
            [1, 1],
            [1, 0],
            [1, 1]
        ]
        papsB = {1: [1, 0], 2: [1, 1], 3: [1, 0], 4: [1, 1]}

        outA = pap2nex(taxa, papsA)
        outB = pap2nex(taxa, papsB)
        
        # nexus format has slightly changed, this touches also the test for the
        # conversion, but the first 100 lines are the important ones, and they
        # are identical with the test string above. we could add teh full
        # current format, but since this may change quickly, it is better to
        # leave the test as this for the moment
        assert nex[:100] == outA[:100] and nex[:100] == outB[:100]

    def test_pap2csv(self):
        csv = """ID	a	b
1	1	0
2	1	1
"""
        paps = {1: [1, 0], 2: [1, 1]}
        taxa = ['a', 'b']
        assert csv == pap2csv(taxa, paps)


class Test_is_constant(TestCase):
    def test_all_present(self):
        assert _is_constant([[2], [1], [4], [5], [3]])
    
    def test_all_absent(self):
        assert _is_constant([0, 0, 0, 0])
    
    def test_not_constant(self):
        assert not _is_constant([0, 0, 0, 0, [8]])
        assert not _is_constant([0, 0, [8], 0, 0])
        assert not _is_constant([[7], [6], 0, 0, 0])
        assert not _is_constant([0, 0, [9], [10], 0])
        assert not _is_constant([[1], [2], [8], 0, 0])
        assert not _is_constant([[1], [2], [8], [3], 0])


class TestWriteNexus(WithTempDir):
    """Tests for `write_nexus`"""
    def setUp(self):
        WithTempDir.setUp(self)
        self.wordlist = Wordlist(test_data('GER.tsv'))
    
    def assertRegexWorkaround(self, a, b):
        # Python has different implementations of assertRegex (>3.1)
        # vs. assertRegexpMatches (2.7, 3.0). This function provides a
        # workaround.
        if sys.version_info >= (3, 2):
            return self.assertRegex(a, b)
        else:
            return self.assertRegexpMatches(a, b)
    
    def test_error_on_unknown_mode(self):
        with self.assertRaises(ValueError):
            write_nexus(self.wordlist, mode='xx')

    def test_error_on_unknown_ref(self):
        with self.assertRaises(KeyError):
            write_nexus(self.wordlist, mode='mrbayes', ref='magic')
    
    def test_mrbayes(self):
        # Use missing="X" parameter to avoid \? in the assertRegex calls below
        nex = write_nexus(self.wordlist, mode='MRBAYES', missing="X",
                filename=text_type(self.tmp_path('test')))
        self.assertIn("NTAX=5 NCHAR=7", nex)
        # mrbayes should have datatype=restriction
        self.assertIn("DATATYPE=RESTRICTION", nex)
        
        # check charblock:
        self.assertRegexWorkaround(nex, r"charset I = 1\-1;")
        self.assertRegexWorkaround(nex, r"charset all = 2\-4;")
        self.assertRegexWorkaround(nex, r"charset ash = 5\-7;")
        
        # check data:
        self.assertRegexWorkaround(nex, r"German\s+1100100")
        self.assertRegexWorkaround(nex, r"English\s+1100XXX")
        self.assertRegexWorkaround(nex, r"Swedish\s+1010010")
        self.assertRegexWorkaround(nex, r"Icelandic\s+1001XXX")
        self.assertRegexWorkaround(nex, r"Norwegian\s+1001001")
      
    def test_splitstree(self):
        # Use missing="X" parameter to avoid \? in the assertRegex calls below
        nex = write_nexus(self.wordlist, mode='SPLITSTREE', missing="X",
                filename=text_type(self.tmp_path('test')))
        self.assertIn("NTAX=5 NCHAR=7", nex)
        # splitstree should have datatype=standard
        self.assertIn("DATATYPE=STANDARD", nex)
        # NO charblock
        assert 'charset' not in nex
        assert 'ASSUMPTIONS' not in nex
        # NO symbols
        assert 'SYMBOLS' not in nex
        # check data:
        self.assertRegexWorkaround(nex, r"German\s+1100100")
        self.assertRegexWorkaround(nex, r"English\s+1100XXX")
        self.assertRegexWorkaround(nex, r"Swedish\s+1010010")
        self.assertRegexWorkaround(nex, r"Icelandic\s+1001XXX")
        self.assertRegexWorkaround(nex, r"Norwegian\s+1001001")
        
    def test_beast(self):
        # Use missing="X" parameter to avoid \? in the assertRegex calls below
        nex = write_nexus(self.wordlist, mode='BEAST', missing="X",
                filename=text_type(self.tmp_path('test')))
        # added one character for ascertainment
        self.assertIn("NTAX=5 NCHAR=8", nex)
        # mrbayes should have datatype=standard
        self.assertIn("DATATYPE=STANDARD", nex)
        
        # check charblock:
        self.assertRegexWorkaround(nex, r"1 _ascertainment,")
        self.assertRegexWorkaround(nex, r"2 I,")
        self.assertRegexWorkaround(nex, r"3 all,")
        self.assertRegexWorkaround(nex, r"4 all,")
        self.assertRegexWorkaround(nex, r"5 all,")
        self.assertRegexWorkaround(nex, r"6 ash,")
        self.assertRegexWorkaround(nex, r"7 ash,")
        self.assertRegexWorkaround(nex, r"8 ash")
        
        # check data:
        self.assertRegexWorkaround(nex, r"German\s+01100100")
        self.assertRegexWorkaround(nex, r"English\s+01100XXX")
        self.assertRegexWorkaround(nex, r"Swedish\s+01010010")
        self.assertRegexWorkaround(nex, r"Icelandic\s+01001XXX")
        self.assertRegexWorkaround(nex, r"Norwegian\s+01001001")
        
    def test_beastwords(self):
        # Use missing="X" parameter to avoid \? in the assertRegex calls below
        nex = write_nexus(self.wordlist, mode='BEASTWORDS', missing="X",
                filename=text_type(self.tmp_path('test')))
        # added three characters for ascertainment
        self.assertIn("NTAX=5 NCHAR=10", nex)
        # mrbayes should have datatype=standard
        self.assertIn("DATATYPE=STANDARD", nex)
        
        # check charblock:
        self.assertRegexWorkaround(nex, r"1 I_ascertainment,")
        self.assertRegexWorkaround(nex, r"2 I,")
        self.assertRegexWorkaround(nex, r"3 all_ascertainment,")
        self.assertRegexWorkaround(nex, r"4 all,")
        self.assertRegexWorkaround(nex, r"5 all,")
        self.assertRegexWorkaround(nex, r"6 all,")
        self.assertRegexWorkaround(nex, r"7 ash_ascertainment,")
        self.assertRegexWorkaround(nex, r"8 ash,")
        self.assertRegexWorkaround(nex, r"9 ash,")
        self.assertRegexWorkaround(nex, r"10 ash")
        
        # check data:
        self.assertRegexWorkaround(nex, r"German\s+0101000100")
        self.assertRegexWorkaround(nex, r"English\s+010100XXXX")
        self.assertRegexWorkaround(nex, r"Swedish\s+0100100010")
        self.assertRegexWorkaround(nex, r"Icelandic\s+010001XXXX")
        self.assertRegexWorkaround(nex, r"Norwegian\s+0100010001")

        # assumptions block
        self.assertRegexWorkaround(nex, r"charset I = 1\-2;")
        self.assertRegexWorkaround(nex, r"charset all = 3\-6;")
        self.assertRegexWorkaround(nex, r"charset ash = 7\-10;")

    def test_traitlab(self):
        nex = write_nexus(self.wordlist, mode='traitlab', missing="X",
                filename=text_type(self.tmp_path('test')))
        # we should lose the FIRST character
        self.assertIn("NTAX=5 NCHAR=6", nex)
        # splitstree should have datatype=standard
        self.assertIn("DATATYPE=STANDARD", nex)
        # NO charblock
        assert 'charset' not in nex
        assert 'ASSUMPTIONS' not in nex
        # NO symbols
        assert 'SYMBOLS' not in nex
        # check data:
        self.assertRegexWorkaround(nex, r"German\s+100100")
        self.assertRegexWorkaround(nex, r"English\s+100XXX")
        self.assertRegexWorkaround(nex, r"Swedish\s+010010")
        self.assertRegexWorkaround(nex, r"Icelandic\s+001XXX")
        self.assertRegexWorkaround(nex, r"Norwegian\s+001001")

    def test_merge_custom_statements(self):
        # this tests for the bug in https://github.com/lingpy/lingpy/issues/340
        import re
        nex = write_nexus(self.wordlist, mode='mrbayes', commands=['test'])
        if len(re.findall(r"BEGIN MRBAYES;", nex, flags=re.IGNORECASE)) == 2:
            raise AssertionError('Duplicate mrbayes block found')

        self.assertRegexWorkaround(nex, r"charset I = 1\-1;")
        self.assertRegexWorkaround(nex, r"charset all = 2\-4;")
        self.assertRegexWorkaround(nex, r"charset ash = 5\-7;")
        self.assertRegexWorkaround(nex, r"test")
        
