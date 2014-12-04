# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-12 12:53
# modified : 2013-11-12 12:53
"""
Test the SCA module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-12"


from six import text_type

from lingpy import Alignments, MSA
import lingpy as lp
from lingpy.tests.util import test_data, WithTempDir


class TestMSA(WithTempDir):
    def test_output(self):
        msa = MSA(test_data('harry.msa'))
        msa.merge = {}  # well. it is a list, but the code apparently wants a dict ...
        fname = text_type(self.tmp_path('test'))
        for fmt in 'msa psa msq html tex'.split():
            msa.output(fileformat=fmt, filename=fname)


class TestAlignments(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.alm = Alignments(test_data('KSL2.qlc'), loans=False)
    
    def test_ipa2tokens(self):
        # iterate over the keys
        for key in self.alm: #.get_list(language="Turkish",flat=True):
            ipa = self.alm[key, 'ipa']
            tokensA = self.alm[key, 'tokensa'].split(' ')
            tokensB = self.alm[key, 'tokensb'].split(' ')

            new_tokensA = lp.ipa2tokens(ipa, merge_vowels=True)
            new_tokensB = lp.ipa2tokens(ipa, merge_vowels=False)
            assert tokensA == new_tokensA
            assert tokensB == new_tokensB

    def test_align(self):
        # align all sequences using standard params
        self.alm.align()

        # iterate and align using the multiple function
        for key,value in self.alm.msa['cogid'].items():
            # first compare simple alignments
            msaA = lp.SCA(value)
            msaB = lp.Multiple(value['seqs'])
            msaB.prog_align()
            assert msaA == msaB

            # now compare with different flag
            msaA = lp.Multiple([self.alm[idx,'tokensb'] for idx in value['ID']])
            msaB = lp.Multiple([''.join(s) for s in value['seqs']],merge_vowels=False)
            msaA.lib_align()
            msaB.lib_align()
            assert msaA == msaB

    def test_get_consensus(self):
        # align all sequences using standard params
        self.alm.align()
        self.alm.get_consensus(consensus="consensus")

        # check whether Turkish strings are identical
        assert self.alm.get_list(
                    language="Turkish",
                    entry="consensus",
                    flat=True
                    ) == \
                            [''.join(x) for x in self.alm.get_list(
                                language="Turkish",
                                entry="tokens",
                                flat=True
                                )
                                ]

    def test_output(self):
        self.alm.align()
        self.alm.output('qlc', filename=text_type(self.tmp_path('test')))
        self.alm.output('html', filename=text_type(self.tmp_path('test')))
