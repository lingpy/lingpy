# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-23 11:00
# modified : 2013-11-23 11:00
"""
Test csv-module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-23"

from lingpy.read.csv import read_asjp
from lingpy.compare.lexstat import LexStat
from lingpy.tests.util import test_data


def test_read_asjp():
    lex = LexStat(read_asjp(
        test_data('asjp_test_list.csv'), family="CELTIC", classification="wls_gen"))
    assert len(lex) == 249
