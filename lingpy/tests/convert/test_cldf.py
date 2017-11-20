from clldutils.testing import WithTempDir
from lingpy.tests.util import test_data

def test_from_cldf():
    from lingpy.basic.wordlist import from_cldf
    wl = from_cldf(test_data('cldf/test-metadata.json'))
    assert wl.width == 29
    assert wl.height == 1
