"""
Test conversions involving strings.
"""
import re

import pytest

from lingpy import rc
from lingpy.algorithm import squareform
from lingpy.align.sca import MSA, get_consensus
from lingpy.basic.wordlist import Wordlist
from lingpy.convert.strings import scorer2str, msa2str, matrix2dst, pap2nex, \
    pap2csv, write_nexus, _is_constant
from lingpy.read import qlc
from lingpy.util import read_text_file


def test_scorer2str(test_data):
    assert scorer2str(rc('dolgo').scorer) == read_text_file(str(test_data / 'dolgo.scorer'))


def test_msa2str(test_data):
    aranger = '{body}{meta}'

    # read msa traditionally into an object
    msa_a = MSA(str(test_data / 'harry.msa'))

    # read msa from dictionary
    msa_b = qlc.read_msa(str(test_data / 'harry.msa'))

    # read msa with IDs
    msa_c = qlc.read_msa(str(test_data / 'harry_with_ids.msa'), ids=True, header=False)

    # we adjust the dataset and the seq_id since otherwise we won't have
    # similar output
    msa_c['seq_id'] = 'test'
    msa_c['dataset'] = 'file'

    # when converting these different objects to string with the same body
    # and the like, they should be identical, so we check this here
    str_a = msa2str(msa_a, _arange=aranger)
    str_b = msa2str(msa_b, _arange=aranger)
    str_c = msa2str(msa_c, _arange=aranger, wordlist=False)

    assert str_a == str_b == str_c

    # we next test for converting with the merging attribute
    str_d = msa2str(msa_c, _arange=aranger, wordlist=True, merge=True)
    str_e = msa2str(msa_c, _arange=aranger, wordlist=True, merge=False)

    # remove tabstops for checking similar strings
    str_d_st = str_d.replace('\t', '')
    str_e_st = str_e.replace('\t', '')

    # get index until 'COLUMN'
    idx = str_d_st.index('COLUMNID')
    assert str_d != str_e and str_d_st[:idx] == str_e_st[:idx]

    # add a consensus string to all msa objects
    consensus_a = get_consensus(MSA(msa_b), gaps=True)
    consensus_b = get_consensus(MSA(msa_c), gaps=True)

    msa_b['consensus'] = consensus_a
    msa_c['consensus'] = consensus_b

    assert msa2str(msa_b) == msa2str(msa_c, wordlist=False)


def test_matrix2dst():
    matrix = squareform([0.5, 0.75, 0.8])

    # we choose same format for taxa as default
    taxa = ['t_1', 't_2', 't_3']

    phyl_a = matrix2dst(matrix, taxa=taxa)
    phyl_b = matrix2dst(matrix)

    assert phyl_a == phyl_b

    phyl_c = matrix2dst(matrix, taxa=taxa, stamp='# Written with joy.')
    phyl_d = matrix2dst(matrix, stamp='# Written with joy.')

    assert phyl_c == phyl_d

    phyl_e = matrix2dst(matrix, taxa=taxa, taxlen=20)
    phyl_f = matrix2dst(matrix, taxlen=30)

    assert 18 * ' ' in phyl_e and 28 * ' ' in phyl_f

    # check for tab-stop output when taxlen is set to 0
    assert matrix2dst(matrix, taxlen=0).count('\t') == 9


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
    paps_a = [
        [1, 0],
        [1, 1],
        [1, 0],
        [1, 1]
    ]
    paps_b = {1: [1, 0], 2: [1, 1], 3: [1, 0], 4: [1, 1]}

    out_a = pap2nex(taxa, paps_a)
    out_b = pap2nex(taxa, paps_b)

    # nexus format has slightly changed, this touches also the test for the
    # conversion, but the first 100 lines are the important ones, and they
    # are identical with the test string above. we could add teh full
    # current format, but since this may change quickly, it is better to
    # leave the test as this for the moment
    assert nex[:100] == out_a[:100] and nex[:100] == out_b[:100]


def test_pap2csv():
    csv = """ID	a	b
1	1	0
2	1	1
"""
    paps = {1: [1, 0], 2: [1, 1]}
    taxa = ['a', 'b']
    assert csv == pap2csv(taxa, paps)


def test_all_present():
    assert _is_constant([[2], [1], [4], [5], [3]])


def test_all_absent():
    assert _is_constant([0, 0, 0, 0])


def test_not_constant():
    assert not _is_constant([0, 0, 0, 0, [8]])
    assert not _is_constant([0, 0, [8], 0, 0])
    assert not _is_constant([[7], [6], 0, 0, 0])
    assert not _is_constant([0, 0, [9], [10], 0])
    assert not _is_constant([[1], [2], [8], 0, 0])
    assert not _is_constant([[1], [2], [8], [3], 0])


@pytest.fixture
def wordlist(test_data):
    return Wordlist(str(test_data / 'GER.tsv'))


def assertRegexWorkaround(string, regex):
    assert re.search(regex, string)


def test_error_on_unknown_mode(wordlist):
    with pytest.raises(ValueError):
        write_nexus(wordlist, mode='xx')


def test_error_on_unknown_ref(wordlist):
    with pytest.raises(KeyError):
        write_nexus(wordlist, mode='mrbayes', ref='magic')


@pytest.fixture
def nexus_factory(wordlist, tmppath):
    def f(**kw):
        # Use missing="X" parameter to avoid \? in the assertRegex calls below
        return write_nexus(wordlist, missing="X", filename=str(tmppath / 'test'), **kw)
    return f


def test_mrbayes(nexus_factory):
    nex = nexus_factory(mode='MRBAYES')
    assert "NTAX=5 NCHAR=7" in nex

    # mrbayes should have datatype=restriction
    assert "DATATYPE=RESTRICTION" in nex

    # check charblock:
    assertRegexWorkaround(nex, r"charset I = 1\-1;")
    assertRegexWorkaround(nex, r"charset all = 2\-4;")
    assertRegexWorkaround(nex, r"charset ash = 5\-7;")

    # check data:
    assertRegexWorkaround(nex, r"German\s+1100100")
    assertRegexWorkaround(nex, r"English\s+1100XXX")
    assertRegexWorkaround(nex, r"Swedish\s+1010010")
    assertRegexWorkaround(nex, r"Icelandic\s+1001XXX")
    assertRegexWorkaround(nex, r"Norwegian\s+1001001")


def test_splitstree(nexus_factory):
    nex = nexus_factory(mode='SPLITSTREE')

    assert "NTAX=5 NCHAR=7" in nex

    # splitstree should have datatype=standard
    assert "DATATYPE=STANDARD" in nex

    # NO charblock
    assert 'charset' not in nex
    assert 'ASSUMPTIONS' not in nex
    # NO symbols
    assert 'SYMBOLS' not in nex
    # check data:
    assertRegexWorkaround(nex, r"German\s+1100100")
    assertRegexWorkaround(nex, r"English\s+1100XXX")
    assertRegexWorkaround(nex, r"Swedish\s+1010010")
    assertRegexWorkaround(nex, r"Icelandic\s+1001XXX")
    assertRegexWorkaround(nex, r"Norwegian\s+1001001")


def test_beast(nexus_factory):
    nex = nexus_factory(mode='BEAST')

    # added one character for ascertainment
    assert "NTAX=5 NCHAR=8" in nex

    # mrbayes should have datatype=standard
    assert "DATATYPE=STANDARD" in nex

    # check charblock:
    assertRegexWorkaround(nex, r"1 _ascertainment,")
    assertRegexWorkaround(nex, r"2 I,")
    assertRegexWorkaround(nex, r"3 all,")
    assertRegexWorkaround(nex, r"4 all,")
    assertRegexWorkaround(nex, r"5 all,")
    assertRegexWorkaround(nex, r"6 ash,")
    assertRegexWorkaround(nex, r"7 ash,")
    assertRegexWorkaround(nex, r"8 ash")

    # check data:
    assertRegexWorkaround(nex, r"German\s+01100100")
    assertRegexWorkaround(nex, r"English\s+01100XXX")
    assertRegexWorkaround(nex, r"Swedish\s+01010010")
    assertRegexWorkaround(nex, r"Icelandic\s+01001XXX")
    assertRegexWorkaround(nex, r"Norwegian\s+01001001")


def test_beastwords(nexus_factory):
    nex = nexus_factory(mode='BEASTWORDS')

    # added three characters for ascertainment
    assert "NTAX=5 NCHAR=10" in nex

    # mrbayes should have datatype=standard
    assert "DATATYPE=STANDARD" in nex

    # check charblock:
    assertRegexWorkaround(nex, r"1 I_ascertainment,")
    assertRegexWorkaround(nex, r"2 I,")
    assertRegexWorkaround(nex, r"3 all_ascertainment,")
    assertRegexWorkaround(nex, r"4 all,")
    assertRegexWorkaround(nex, r"5 all,")
    assertRegexWorkaround(nex, r"6 all,")
    assertRegexWorkaround(nex, r"7 ash_ascertainment,")
    assertRegexWorkaround(nex, r"8 ash,")
    assertRegexWorkaround(nex, r"9 ash,")
    assertRegexWorkaround(nex, r"10 ash")

    # check data:
    assertRegexWorkaround(nex, r"German\s+0101000100")
    assertRegexWorkaround(nex, r"English\s+010100XXXX")
    assertRegexWorkaround(nex, r"Swedish\s+0100100010")
    assertRegexWorkaround(nex, r"Icelandic\s+010001XXXX")
    assertRegexWorkaround(nex, r"Norwegian\s+0100010001")

    # assumptions block
    assertRegexWorkaround(nex, r"charset I = 1\-2;")
    assertRegexWorkaround(nex, r"charset all = 3\-6;")
    assertRegexWorkaround(nex, r"charset ash = 7\-10;")


def test_traitlab(nexus_factory):
    nex = nexus_factory(mode='traitlab')

    # we should lose the FIRST character
    assert "NTAX=5 NCHAR=6" in nex

    # splitstree should have datatype=standard
    assert "DATATYPE=STANDARD" in nex

    # NO charblock
    assert 'charset' not in nex
    assert 'ASSUMPTIONS' not in nex
    # NO symbols
    assert 'SYMBOLS' not in nex
    # check data:
    assertRegexWorkaround(nex, r"German\s+100100")
    assertRegexWorkaround(nex, r"English\s+100XXX")
    assertRegexWorkaround(nex, r"Swedish\s+010010")
    assertRegexWorkaround(nex, r"Icelandic\s+001XXX")
    assertRegexWorkaround(nex, r"Norwegian\s+001001")


def test_merge_custom_statements(nexus_factory):
    # this tests for the bug in https://github.com/lingpy/lingpy/issues/340
    nex = nexus_factory(mode='mrbayes', commands=['test'])
    assert len(re.findall(r"BEGIN MRBAYES;", nex, flags=re.IGNORECASE)) < 2

    assertRegexWorkaround(nex, r"charset I = 1\-1;")
    assertRegexWorkaround(nex, r"charset all = 2\-4;")
    assertRegexWorkaround(nex, r"charset ash = 5\-7;")
    assertRegexWorkaround(nex, r"test")
