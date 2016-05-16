from lingpy.read.qlc import read_msa, reduce_alignment, read_qlc
from lingpy.align import MSA
from lingpy.tests.util import test_data


def test_read_msa():
    msa = MSA(read_msa(test_data('harry.msa')))
    assert hasattr(msa, 'seqs')

def test_normalize_alignment():

    msa = MSA(read_msa(test_data('harry_unnormal.msa')))

    for line in msa.alignment[1:]:
        assert len(line) == len(msa.alignment[0])

def test_reduce_msa():

    msa = MSA(read_msa(test_data('test_reduce.msa')))
    reduced_alignment = reduce_alignment(msa.alignment)
    for i,line in enumerate(reduced_alignment):
        assert len(line) == 4 and \
                ''.join(line) == ''.join(
                        msa.alignment[i])[:msa.alignment[i].index('(')]

def test_read_qlc():

    qlc = read_qlc(test_data('read_qlc.qlc'))
