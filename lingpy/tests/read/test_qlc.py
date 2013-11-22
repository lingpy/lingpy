# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-13 13:19
# modified : 2013-11-13 13:19
"""
Test qlc parsing module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-13"

import os
from lingpy.read.qlc import read_msa, _list2msa
from lingpy.align import MSA
from lingpy.settings import rcParams

def test_read_msa():

    msa_dict = read_msa(os.path.join(
        rcParams['_path'],
        'tests',
        'test_data',
        'harry.msa')
        )

    msa = MSA(msa_dict)
    assert hasattr(msa,'seqs') == True

