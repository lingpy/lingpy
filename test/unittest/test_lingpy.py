# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-21 12:36
# modified : 2013-07-21 12:36
"""
Carry out test of the alignment module.
"""

__author__="Johann-Mattis List"
__date__="2013-07-21"

# import lingpy library
import lingpy as lp
import os

# import unittest
import unittest

class TestAlignments(unittest.TestCase):

    def setUp(self):

        self.alm = lp.Alignments(
                os.path.join('data','kessler.qlc'),
                loans = False
                )
    
    def test_ipa2tokens(self):
        
        # iterate over the keys
        for key in self.alm: #.get_list(language="Turkish",flat=True):

            ipa = self.alm[key,'ipa']
            tokensA = self.alm[key,'tokensa'].split(' ')
            tokensB = self.alm[key,'tokensb'].split(' ')

            new_tokensA = lp.ipa2tokens(ipa,merge_vowels=True)
            new_tokensB = lp.ipa2tokens(ipa,merge_vowels=False)
            self.assertEqual(tokensA,new_tokensA)
            self.assertEqual(tokensB,new_tokensB)

    def test_align(self):
        
        # align all sequences using standard params
        self.alm.align()

        # iterate and align using the multiple function
        for key,value in self.alm.msa['cogid'].items():

            # first compare simple alignments
            msaA = lp.SCA(value)
            msaB = lp.Multiple(value['seqs'])
            msaB.prog_align()
            self.assertEqual(msaA,msaB)

            # now compare with different flag
            msaA = lp.Multiple([self.alm[idx,'tokensb'] for idx in value['ID']])
            msaB = lp.Multiple([''.join(s) for s in value['seqs']],merge_vowels=False)
            msaA.lib_align()
            msaB.lib_align()
            self.assertEqual(msaA,msaB)


    def test_get_consensus(self):

        
        # align all sequences using standard params
        self.alm.align()
        self.alm.get_consensus(consensus="consensus")

        # check whether Turkish strings are identical
        self.assertEqual(
                self.alm.get_list(
                    language="Turkish",
                    entry="consensus",
                    flat=True
                    ),
                [''.join(x) for x in self.alm.get_list(
                    language="Turkish",
                    entry="tokens",
                    flat=True
                    )
                    ]
                )
        


# )


class TestLexStat(unittest.TestCase):
    
    def setUp(self):

        self.lex = lp.LexStat(
                os.path.join(
                    "data",
                    "kessler.qlc"
                    )
                )

    def test_get_scorer(self):
            
        self.lex.get_scorer()
        self.assertEqual(
                hasattr(self.lex,"cscorer"),
                True
                )

    def test_cluster(self):
        
        self.lex.get_scorer()
        self.lex.cluster(method="lexstat",threshold=0.7)
        self.lex.cluster(method="edit-dist",threshold=0.7)
        self.lex.cluster(method="turchin",threshold=0.7)

        self.assertEqual(
                ('scaid' in self.lex.header and
                    'lexstatid' in self.lex.header and
                    'editid' in self.lex.header and
                    'turchinid' in self.lex.header),
                True
                )


if __name__ == "__main__":
    unittest.main()

