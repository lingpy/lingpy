# *-* coding: utf-8 *-* 
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from ..util import test_data
from nose.tools import assert_raises
from six import text_type
from lingpy.basic.wordlist import Wordlist
from lingpy.sequence.profile import simple_profile, context_profile


def test_simple_profile():
    wl = Wordlist(test_data('KSL.qlc'))
    prf = list(simple_profile(wl))
    assert ('a', 'a', '507', 'U+0061') in prf

def test_context_profile():
    wl = Wordlist(test_data('KSL.qlc'))
    prf = list(context_profile(wl))
    assert '376' in prf[0] # first line of profile

