# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-13 13:19
# modified : 2013-11-13 13:19
"""
Test qlc parsing module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-13"

from lingpy.read.qlc import read_msa
from lingpy.align import MSA
from lingpy.tests.util import test_data


def test_read_msa():
    msa = MSA(read_msa(test_data('harry.msa')))
    assert hasattr(msa, 'seqs')
