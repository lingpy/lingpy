import pytest

from lingpy.align.sca import MSA, PSA
from lingpy.evaluate.apa import EvalPSA, EvalMSA


def test_EvalPSA(test_data, tmp_path):
    obj = EvalPSA(
        PSA(str(test_data / 'harry_potter.psa')),
        PSA(str(test_data / 'harry_potter_misaligned.psa')))
    obj.c_score()
    obj.r_score()
    obj.sp_score()
    obj.jc_score()
    obj.diff(filename=str(tmp_path / 'test_EvalPSA.diff'))


def test_EvalMSA(test_data):
    msa = MSA(str(test_data / 'harry.msa'))
    msa2 = MSA(str(test_data / 'harryp.msa'))

    for test in [msa, msa2]:
        obj = EvalMSA(msa, test)
        for mode in range(1, 5):
            obj.c_score(mode=mode)
            if hasattr(obj, 'pic'):
                del obj.pic  # pragma: no cover
        with pytest.raises(ValueError):
            obj.c_score(10)
        res = obj.r_score()
        if test == msa:
            assert res == pytest.approx(1.0)
        obj.sp_score()
        obj.jc_score()
        obj.check_swaps()
