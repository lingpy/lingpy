# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-11-23 11:00
# modified : 2013-11-23 11:00
"""
Test csv-module.
"""

__author__="Johann-Mattis List"
__date__="2013-11-23"

from lingpy import rc
from lingpy.read.csv import *
from lingpy.compare.lexstat import LexStat

def test_read_asjp():

    p = rc('test_path')+'asjp_test_list.csv'
    
    data = read_asjp(p, family="CELTIC", classification="wls_gen")

    lex = LexStat(data)
    assert len(lex) == 249
    
