# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-12-09 09:40
# modified : 2014-12-09 09:40
"""
Test conversions involving strings.
"""

__author__="Johann-Mattis List"
__date__="2014-12-09"

import lingpy
from lingpy.convert.strings import *
from ..util import *
import os

def test_scorer2str():
    """
    Test conversion of scorers to strings.
    """
    
    # get scorer from dolgo (not so many chars, easier to test)
    scorerA = scorer2str(lingpy.rc('dolgo').scorer)

    # get scoring matrix for dolgo
    scorerB = lingpy.util.read_text_file(os.path.join(test_data(),'dolgo.scorer'))

    assert scorerA == scorerB

def test_msa2str():

    aranger = '{body}{meta}'

    
    # read msa traditionally into an object
    msaA = lingpy.MSA(os.path.join(test_data(),'harry.msa'))
    
    # read msa from dictionary
    msaB = lingpy.read.qlc.read_msa(os.path.join(test_data(),'harry.msa'))

    # read msa with IDs
    msaC = lingpy.read.qlc.read_msa(os.path.join(test_data(),'harry_with_ids.msa'),
            ids=True, header=False)

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
    strDst = strD.replace('\t','')
    strEst = strE.replace('\t','')

    # get index until 'COLUMN'
    idx = strDst.index('COLUMNID')
    assert strD != strE and strDst[:idx] == strEst[:idx]

    # add a consensus string to all msa objects
    consensusA = lingpy.align.sca.get_consensus(lingpy.align.sca.MSA(msaB), gaps=True)
    consensusB = lingpy.align.sca.get_consensus(lingpy.align.sca.MSA(msaC), gaps=True)

    msaB['consensus'] = consensusA
    msaC['consensus'] = consensusB

    assert msa2str(msaB) == msa2str(msaC, wordlist=False)

def test_matrix2dst():

    matrix = lingpy.algorithm.squareform([0.5, 0.75, 0.8])
    
    # we choose same format for taxa as default
    taxa = ['t_1', 't_2', 't_3']

    phylA = matrix2dst(matrix, taxa=taxa)
    phylB = matrix2dst(matrix)

    assert phylA == phylB
    
    phylC = matrix2dst(matrix, taxa=taxa, stamp='# Written with joy.')
    phylD = matrix2dst(matrix, stamp = '# Written with joy.')

    assert phylC == phylD

    phylE = matrix2dst(matrix, taxa=taxa, taxlen=20)
    phylF = matrix2dst(matrix, taxlen=30)

    assert 18 * ' ' in phylE and 28 * ' ' in phylF
    
    # check for tab-stop output when taxlen is set to 0

    phylG = matrix2dst(matrix, taxlen=0)
    assert phylG.count('\t') == 9


def test_pap2nex():

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
            [1,0],
            [1,1],
            [1,0],
            [1,1]
            ]
    papsB = {1: [1,0], 2: [1,1], 3:[1,0], 4:[1,1]}

    outA = pap2nex(taxa, papsA)
    outB = pap2nex(taxa, papsB)
    
    assert nex == outA and nex == outB

def test_pap2csv():
    
    csv = """ID	a	b
1	1	0
2	1	1
"""

    paps = {1: [1,0], 2: [1,1]}

    taxa = ['a','b']

    assert csv == pap2csv(taxa, paps)


