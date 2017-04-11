from __future__ import (
        unicode_literals, print_function, absolute_import, division)

from lingpy.tests.util import test_data
from lingpy.basic.wordlist import Wordlist
from lingpy.compare.sanity import mutual_coverage

def test_mutual_coverage():
    wl = Wordlist(test_data('KSL5.qlc'))
    mc = mutual_coverage(wl, 0)
    assert len(mc[2][0]) == wl.width
    mc = mutual_coverage(wl, 3)
    assert len(mc[3][0]) == 3
    assert not mutual_coverage(wl, 4)
