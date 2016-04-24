from lingpy.convert.html import *
from lingpy import *
from lingpy.tests.util import WithTempDir
from lingpy.tests.util import test_data
from six import text_type

class Tests(WithTempDir):


    def test_color_range(self):

        crang = colorRange(10)
        assert sorted(set([x[0] for x in crang])) == ['#']
        assert len(set(crang)) == len(crang)
    
    def test_alm2html(self):

        alm2html(test_data('alm.alm'), filename=text_type(self.tmp_path('alm')))

    def test_msa2html(self):

        msa = MSA(test_data('harry.msq'))
        msa.prog_align()
        msa.output('html', filename=text_type(self.tmp_path('alm')))

    def test_strings_and_tokens2html(self):

        tokens2html(list('haXy'))
        tokens2html(list('hary'), swaps=[1,2])
        tokens2html(list('hary'), tax_len=20)
        string2html('English', list('hary'))
        string2html('English', list('haXy'), swaps=[1,2])
        string2html('English', list('hary'), tax_len=20)

    def test_psa2html(self):

        psa2html(test_data('harry_potter.psa'),
                filename=text_type(self.tmp_path('alm')))
        psa2html(test_data('harry_potter_bad.psa'),
                filename=text_type(self.tmp_path('alm')))
    
