from lingpy.basic.wordlist import Wordlist
from lingpy.sequence.profile import simple_profile, context_profile


def test_simple_profile(test_data):
    wl = Wordlist(str(test_data / 'KSL6.qlc'))
    prf = list(simple_profile(wl))
    assert ('a', 'a', '7', 'U+0061') in prf
    prf = list(simple_profile(wl, clts={'a': 'A'}))
    assert prf[0][1] == 'A'


def test_context_profile(test_data):
    wl = Wordlist(str(test_data / 'KSL6.qlc'))
    prf = list(context_profile(wl))
    assert prf[2][-2] == '4'  # first line of profile
    prf = list(context_profile(wl, clts={'a': 'A'}))
    assert prf[2][1] == 'A'
