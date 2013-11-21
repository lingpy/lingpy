# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-12 12:53
# modified : 2013-11-12 12:53
"""
Test the SCA module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-12"

import os
import unittest
from lingpy import Alignments
from lingpy.settings import rcParams
import lingpy as lp

class TestAlignments(object):

    def setup(self):

        self.alm = Alignments(
                os.path.join(
                    rcParams['_path'],
                    'tests',
                    'test_data',
                    'KSL2.qlc'
                    ),
                loans=False
                )
    
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
        
        try:
            self.alm.align()
            self.alm.output('qlc', filename='test')
            self.alm.output('html', filename='test')

            assert True

        except:
            assert False
