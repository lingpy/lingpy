from lingpy.tests.util_testing import WithTempDir
from lingpy.tests.util import test_data

def test_from_cldf():
    from lingpy.basic.wordlist import from_cldf
    wl = from_cldf(test_data('cldf/test-metadata.json'), language='Name',
            concept='Name', concepticon="Concepticon_ID",
            glottocode='glottocode')
    assert wl.width == 29
    assert wl.height == 1
