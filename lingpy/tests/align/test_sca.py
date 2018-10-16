"""
Test the SCA module.
"""
from __future__ import unicode_literals

from itertools import product

from six import text_type

import lingpy as lp
from lingpy import Alignments, MSA, PSA, LexStat
from lingpy.tests.util import test_data
from lingpy.tests.util_testing import WithTempDir
from lingpy.util import write_text_file


class TestPSA(WithTempDir):
    def test_output(self):
        fpsa = self.tmp_path('test.psa')
        write_text_file(fpsa, '\n')
        psa = PSA(text_type(fpsa))
        fname = text_type(self.tmp_path('test'))
        psa.output(fileformat='psa', filename=fname)

        psq = self.tmp_path('test.psq')
        write_text_file(psq, '\n')
        psa = PSA(text_type(psq))
        fname = text_type(self.tmp_path('test'))
        psa.output(fileformat='psq', filename=fname)

        psa = PSA(text_type(test_data('harry_potter.psa')))
        psa.align()
        psa.output(fileformat="psa", filename=fname, scores=True)
        psa.output(fileformat="psq", filename=fname)


class TestMSA(WithTempDir):
    def test_output(self):
        msa = MSA(test_data('harry.msa'))
        msa.ipa2cls()
        # well. it is a list, but the code apparently wants a dict ...
        msa.merge = {'a': 'x', 'b': 'x'}
        fname = text_type(self.tmp_path('test'))
        for fmt in 'msa psa msq html tex'.split():
            for s, u in product([True, False], [True, False]):
                msa.output(fileformat=fmt, filename=fname, sorted_seqs=s,
                           unique_seqs=u)


class TestAlignments(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.alm = Alignments(test_data('KSL2.qlc'), loans=False,
                              _interactive=False)
        self.alm.align()

    def test_ipa2tokens(self):
        # iterate over the keys
        for key in self.alm:  # get_list(language="Turkish",flat=True):
            ipa = self.alm[key, 'ipa']
            tokens_a = self.alm[key, 'tokensa'].split(' ')
            tokens_b = self.alm[key, 'tokensb'].split(' ')

            new_tokens_a = lp.ipa2tokens(ipa, merge_vowels=True,
                                         merge_geminates=False)
            new_tokens_b = lp.ipa2tokens(ipa, merge_vowels=False,
                                         merge_geminates=False)
            assert tokens_a == new_tokens_a
            assert tokens_b == new_tokens_b

    def test_align(self):
        self.alm.add_entries('cugid', self.alm._ref, lambda x: text_type(x))
        self.alm.add_alignments(ref="cugid")

        # align all sequences using standard params
        self.alm.align(ref="cugid", alignment="alignment2")
        assert (self.alm.msa["cugid"]["1"]["ID"] ==
                self.alm.msa["cogid"][1]["ID"])

        # iterate and align using the multiple function
        for key, value in self.alm.msa['cogid'].items():
            # first compare simple alignments
            msa_a = lp.SCA(value)
            msa_b = lp.Multiple(value['seqs'])
            msa_b.prog_align()
            assert msa_a == msa_b

            # now compare with different flag
            msa_a = lp.Multiple([self.alm[idx, 'tokensb']
                                 for idx in value['ID']])
            msa_b = lp.Multiple([''.join(s) for s in value['seqs']],
                                merge_vowels=False)
            msa_a.lib_align()
            msa_b.lib_align()
            assert msa_a == msa_b

    def test_get_consensus(self):
        # align all sequences using standard params
        self.alm.get_consensus(consensus="consensus", classes=True)
        self.alm.get_consensus(consensus="consensus")

        # check whether Turkish strings are identical
        self.assertEqual(
            self.alm.get_list(language="Turkish", entry="consensus", flat=True),
            [''.join(x) for x in
             self.alm.get_list(language="Turkish", entry="tokens", flat=True)])

    def test_get_confidence(self):
        lex = LexStat(test_data('KSL3.qlc'))
        tmp_dict = dict([(k, lex[k, 'numbers']) for k in lex])
        self.alm.add_entries('numbers', tmp_dict, lambda x: x)
        # Run get_confidence to populate the output variable.
        # TODO: Check and document side-effects of this.
        _ = self.alm.get_confidence(lex.rscorer, ref='cogid')
        self.alm.output('html', filename=text_type(self.tmp_path('alm')),
                        confidence=True)

    def test_output(self):
        self.alm.output('tsv', filename=text_type(self.tmp_path('test')))
        self.alm.output('html', filename=text_type(self.tmp_path('test')))


def test_get_consensus():
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

def test_partial_alignments_with_lexstat():
    lex = lp.LexStat(test_data('test-partial-alignments.tsv'), segments='tokens')
    alms = lp.Alignments(test_data('test-partial-alignments.tsv'), fuzzy=True,
            ref='cogids', sonar=True, segments='tokens')
    alms.align(scorer=lex.bscorer)
    assert '-' in alms.msa['cogids'][12]['alignment'][-1]
