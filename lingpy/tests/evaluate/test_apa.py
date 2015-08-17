# coding: utf8
from __future__ import unicode_literals

from lingpy.align.sca import MSA, PSA
from lingpy.tests.util import test_data, WithTempDir


class Tests(WithTempDir):
    def test_EvalPSA(self):
        from lingpy.evaluate.apa import EvalPSA

        obj = EvalPSA(
            PSA(test_data('harry_potter.psa')),
            PSA(test_data('harry_potter_misaligned.psa')))
        obj.c_score()
        obj.r_score()
        obj.sp_score()
        obj.jc_score()
        obj.diff(filename='%s' % self.tmp_path('test_EvalPSA.diff'))

    def test_EvalMSA(self):
        from lingpy.evaluate.apa import EvalMSA

        msa = MSA(test_data('harry.msa'))
        msa2 = MSA(test_data('harryp.msa'))

        for test in [msa, msa2]:
            obj = EvalMSA(msa, test)
            for mode in range(1, 5):
                obj.c_score(mode=mode)
                if hasattr(obj, 'pic'):
                    del obj.pic
            self.assertRaises(ValueError, obj.c_score, 10)
            res = obj.r_score()
            if test == msa:
                self.assertAlmostEquals(res, 1.0)
            obj.sp_score()
            obj.jc_score()
            obj.check_swaps()
