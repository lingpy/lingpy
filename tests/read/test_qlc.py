import pytest

from lingpy.align import MSA
from lingpy.read.qlc import read_msa, reduce_alignment, read_qlc
from lingpy.thirdparty.cogent.tree import TreeNode


@pytest.mark.parametrize(
    'input,output',
    [
        (['(a b c', 'c d) e'], [['-', 'a', ' ', 'b', ' ', 'c'], ['c', ' ', 'd', '-', ' ', 'e']]),
    ]
)
def test_reduce_alignment(input, output):
    assert reduce_alignment(input) == output


def test_read_msa(test_data):
    msa = MSA(read_msa(str(test_data / 'harry.msa')))
    assert hasattr(msa, 'seqs')


def test_normalize_alignment(test_data):
    msa = MSA(read_msa(str(test_data / 'harry_unnormal.msa')))

    for line in msa.alignment[1:]:
        assert len(line) == len(msa.alignment[0])


def test_reduce_msa(test_data):
    msa = MSA(read_msa(str(test_data / 'test_reduce.msa')))
    reduced_alignment = reduce_alignment(msa.alignment)
    for i, line in enumerate(reduced_alignment):
        assert len(line) == 4 and \
                ''.join(line) == ''.join(
                        msa.alignment[i])[:msa.alignment[i].index('(')]


def test_read_qlc(test_data):
    _ = read_qlc(str(test_data / 'read_qlc.qlc'))


def test_read_qlc_complex(tmp_path):
    p = tmp_path / 'test.qlc'
    p.write_text("""\
<json id="x">
{"y": 5}
</json>
<nwk>
(A,B)C;
</nwk>
<csv id="z" dtype="int" ncol="3">
a\t4\t5
</csv>
ID  NAME
x   z
""", encoding='utf8')
    res = read_qlc(str(p))
    assert res['x']['y'] == 5
    assert isinstance(res['trees']['1'], TreeNode)
    assert res['z']['a'] == [4, 5]
