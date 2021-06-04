from lingpy.align.sca import MSA
from lingpy.convert.html import *


def test_color_range():
    crang = colorRange(10)
    assert sorted(set([x[0] for x in crang])) == ['#']
    assert len(set(crang)) == len(crang)


def test_alm2html(test_data, tmp_path):
    alm2html(str(test_data / 'alm.alm'), filename=str(tmp_path / 'alm'))


def test_msa2html(test_data, tmp_path):
    msa = MSA(str(test_data / 'harry.msq'))
    msa.prog_align()
    msa.output('html', filename=str(tmp_path / 'alm'))


def test_strings_and_tokens2html():
    tokens2html(list('haXy'))
    tokens2html(list('hary'), swaps=[1, 2])
    tokens2html(list('hary'), tax_len=20)
    string2html('English', list('hary'))
    string2html('English', list('haXy'), swaps=[1, 2])
    string2html('English', list('hary'), tax_len=20)


def test_psa2html(test_data, tmp_path):
    psa2html(str(test_data / 'harry_potter.psa'), filename=str(tmp_path / 'alm'))
    psa2html(str(test_data / 'harry_potter_bad.psa'), filename=str(tmp_path / 'alm'))
