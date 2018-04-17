"""
Basic tests for the Phylip module.
"""

from lingpy.read.phylip import *
from lingpy.tests.util import test_data

from unittest import TestCase

MATRIX = """4
German	0.0	0.8	0.4	0.7
English	0.8	0.0	0.5	0.8
French	0.4	0.5	0.0	0.3
Norwegian	0.7	0.8	0.3	0.0"""

DOLGO = """+	0.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-100.00	-5.00	-100.00
0	-100.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00
1	-100.00	0.00	2.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	0.00	1.00
H	-100.00	0.00	-20.00	10.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	-10.00	0.00	0.00	-20.00
J	-100.00	0.00	-20.00	0.00	10.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	-10.00	0.00	0.00	-20.00
K	-100.00	0.00	-20.00	0.00	0.00	10.00	0.00	0.00	0.00	0.00	0.00	0.00	-10.00	0.00	0.00	-20.00
M	-100.00	0.00	-20.00	0.00	0.00	0.00	10.00	0.00	0.00	0.00	0.00	0.00	-10.00	0.00	0.00	-20.00
N	-100.00	0.00	-20.00	0.00	0.00	0.00	0.00	10.00	0.00	0.00	0.00	0.00	-10.00	0.00	0.00	-20.00
P	-100.00	0.00	-20.00	0.00	0.00	0.00	0.00	0.00	10.00	0.00	0.00	0.00	-10.00	0.00	0.00	-20.00
R	-100.00	0.00	-20.00	0.00	0.00	0.00	0.00	0.00	0.00	10.00	0.00	0.00	-10.00	0.00	0.00	-20.00
S	-100.00	0.00	-20.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	10.00	0.00	-10.00	0.00	0.00	-20.00
T	-100.00	0.00	-20.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	10.00	-10.00	0.00	0.00	-20.00
V	-100.00	0.00	-20.00	-10.00	-10.00	-10.00	-10.00	-10.00	-10.00	-10.00	-10.00	-10.00	5.00	-10.00	0.00	-20.00
W	-100.00	0.00	-20.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	-10.00	10.00	0.00	-20.00
X	-5.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00
_	-100.00	0.00	1.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	-20.00	0.00	2.00"""


class Tests(TestCase):
    def setUp(self):
        pass

    def test_read_dst(self):
        t1, m1 = read_dst(test_data('phylip_basic.dst'))
        t2, m2 = read_dst(test_data('phylip_tabstop.dst'), taxlen=0)
        t3, m3 = read_dst(MATRIX, taxlen=0)

        assert t1 == t2 == t3

        ma0 = sum([m[0] for m in m1])  # 1.9
        ma1 = sum([m[1] for m in m1])  # 2.1
        ma2 = sum([m[2] for m in m1])  # 1.2
        ma3 = sum([m[3] for m in m1])  # 1.8
        mb0 = sum([m[0] for m in m2])  # 1.9
        mb1 = sum([m[1] for m in m2])  # 2.1
        mb2 = sum([m[2] for m in m2])  # 1.2
        mb3 = sum([m[3] for m in m2])  # 1.8

        assert round(ma0, 2) == round(mb0, 2) == 1.9
        assert round(ma1, 2) == round(mb1, 2) == 2.1
        assert round(ma2, 2) == round(mb2, 2) == 1.2
        assert round(ma3, 2) == round(mb3, 2) == 1.8

    def test_read_scorer(self):
        scorer = read_scorer(test_data('dolgo.matrix'))

        assert sorted(scorer.chars2int)[0] == '+'
        for letter in 'PTKRSM':
            assert scorer[letter, 'V'] == -10

        assert max(scorer.chars2int.values()) == 15

        # add scorer from string
        scorer2 = read_scorer(DOLGO)
        assert sorted(scorer.chars2int) == sorted(scorer2.chars2int)
