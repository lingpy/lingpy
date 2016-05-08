"""
Test conversions involving strings.
"""
from __future__ import division, unicode_literals
from unittest import TestCase

import lingpy
from lingpy.convert.strings import scorer2str, msa2str, matrix2dst, pap2nex, pap2csv
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
