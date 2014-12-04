# *-* coding: utf-8 *-*
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-12-02 15:37
# modified : 2014-12-02 15:37
"""
Testing multiple module.
"""

__author__="Johann-Mattis List"
__date__="2014-12-02"

import os
import unittest
from lingpy.align import Pairwise, pw_align, nw_align, edit_dist, sw_align,\
        we_align, structalign, turchin


class TestPairwise(object):

    def setup(self):

        self.pair = Pairwise('waldemar','vladimir')

    def test_align(self):

        self.pair.align()

        assert ''.join(self.pair.alignments[0][0]) == 'wal-demar'
        assert ''.join(self.pair.alignments[0][1]) == 'v-ladimir'
    
def test_pw_align():

    assert pw_align('waldemar','vladimir', mode='local')[2] == 1

def test_nw_align():
    
    assert nw_align('waldemar','vladimir')[2] == -1

def test_sw_align():

    assert sw_align('waldemar','wald')[2] == 4

def test_we_align():

    assert we_align('waldemar','marwalde')[1][2] == 3

def test_structalign():
    
    assert len(structalign('waldemar','vladimir')[0]) == 3

def test_turchin():

    assert turchin('waldemar','vladimir') == 0


