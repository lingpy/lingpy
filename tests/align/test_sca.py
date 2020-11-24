"""
Test the SCA module.
"""
from itertools import product

import pytest

import lingpy as lp
from lingpy import Alignments, MSA, PSA, LexStat
from lingpy.util import write_text_file


def test_output(tmppath, test_data):
    fpsa = tmppath / 'test.psa'
    write_text_file(fpsa, '\n')
    psa = PSA(str(fpsa))
    fname = str(tmppath / 'test')
    psa.output(fileformat='psa', filename=fname)

    psq = tmppath / 'test.psq'
    write_text_file(psq, '\n')
    psa = PSA(str(psq))
    fname = str(tmppath / 'test')
    psa.output(fileformat='psq', filename=fname)

    psa = PSA(str(test_data / 'harry_potter.psa'))
    psa.align()
    psa.output(fileformat="psa", filename=fname, scores=True)
    psa.output(fileformat="psq", filename=fname)


def test_output2(test_data, tmppath):
    msa = MSA(str(test_data / 'harry.msa'))
    msa.ipa2cls()
    # well. it is a list, but the code apparently wants a dict ...
    msa.merge = {'a': 'x', 'b': 'x'}
    fname = str(tmppath / 'test')
    for fmt in 'msa psa msq html tex'.split():
        for s, u in product([True, False], [True, False]):
            msa.output(fileformat=fmt, filename=fname, sorted_seqs=s, unique_seqs=u)


@pytest.fixture
def alm(test_data):
    a = Alignments(str(test_data / 'KSL2.qlc'), loans=False,_interactive=False)
    a.align()
    return a


def test_ipa2tokens(alm):
    # iterate over the keys
    for key in alm:  # get_list(language="Turkish",flat=True):
        ipa = alm[key, 'ipa']
        tokens_a = alm[key, 'tokensa'].split(' ')
        tokens_b = alm[key, 'tokensb'].split(' ')

        new_tokens_a = lp.ipa2tokens(ipa, merge_vowels=True, merge_geminates=False)
        new_tokens_b = lp.ipa2tokens(ipa, merge_vowels=False, merge_geminates=False)
        assert tokens_a == new_tokens_a
        assert tokens_b == new_tokens_b


def test_align(alm):
    alm.add_entries('cugid', alm._ref, lambda x: str(x))
    alm.add_alignments(ref="cugid")

    # align all sequences using standard params
    alm.align(ref="cugid", alignment="alignment2")
    assert alm.msa["cugid"]["1"]["ID"] == alm.msa["cogid"][1]["ID"]

    # iterate and align using the multiple function
    for key, value in alm.msa['cogid'].items():
        # first compare simple alignments
        msa_a = lp.SCA(value)
        msa_b = lp.Multiple(value['seqs'])
        msa_b.prog_align()
        assert msa_a == msa_b

        # now compare with different flag
        msa_a = lp.Multiple([alm[idx, 'tokensb'] for idx in value['ID']])
        msa_b = lp.Multiple([''.join(s) for s in value['seqs']], merge_vowels=False)
        msa_a.lib_align()
        msa_b.lib_align()
        assert msa_a == msa_b


def test_get_consensus(alm):
    # align all sequences using standard params
    alm.get_consensus(consensus="consensus", classes=True)
    alm.get_consensus(consensus="consensus")

    # check whether Turkish strings are identical
    assert alm.get_list(language="Turkish", entry="consensus", flat=True) == \
        [''.join(x) for x in
         alm.get_list(language="Turkish", entry="tokens", flat=True)]


def test_get_confidence(test_data, alm, tmppath):
    lex = LexStat(str(test_data / 'KSL3.qlc'))
    tmp_dict = dict([(k, lex[k, 'numbers']) for k in lex])
    alm.add_entries('numbers', tmp_dict, lambda x: x)
    # Run get_confidence to populate the output variable.
    # TODO: Check and document side-effects of this.
    _ = alm.get_confidence(lex.rscorer, ref='cogid')
    alm.output('html', filename=str(tmppath / 'alm'), confidence=True)


def test_output3(alm, tmppath):
    alm.output('tsv', filename=str(tmppath / 'test'))
    alm.output('html', filename=str(tmppath / 'test'))


def test_get_consensus2():
    strings = ['harry', 'harald', 'gari']
    classes = lp.algorithm.misc.transpose([list('H--ARY'), list('HARALT'),
                                           list('KAR--I')])

    msa = lp.align.multiple.mult_align(strings)
    cons = lp.align.sca.get_consensus(msa)
    cons2 = lp.align.sca.get_consensus(msa, gaps=True)
    cons3 = lp.align.sca.get_consensus(msa, classes=classes)
    cons4 = lp.align.sca.get_consensus(msa, local="peaks")
    cons5 = lp.align.sca.get_consensus(msa, local="gaps")

    assert cons == [x for x in cons2 if x != '-']
    assert cons3[:2] == cons[:2]
    assert cons4[0] == 'h'
    assert cons5[0] == 'h'


def test_partial_alignments_with_lexstat(test_data):
    lex = lp.LexStat(str(test_data / 'test-partial-alignments.tsv'), segments='tokens')
    alms = lp.Alignments(str(test_data / 'test-partial-alignments.tsv'), fuzzy=True,
            ref='cogids', sonar=True, segments='tokens')
    alms.align(scorer=lex.bscorer)
    assert '-' in alms.msa['cogids'][12]['alignment'][-1]
