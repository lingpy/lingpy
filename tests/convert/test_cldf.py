from lingpy.basic.wordlist import from_cldf


def test_from_cldf(test_data):
    wl = from_cldf(str(test_data / 'cldf/test-metadata.json'), language='Name',
                   concept='Name', concepticon="Concepticon_ID",
                   glottocode='glottocode')

    assert wl.width == 29
    assert wl.height == 1
    assert wl.entries[0] == 'alignment'
    assert wl.cols[0] == 'Anuta'
    assert wl.cols[28] == 'Vaeakau-Taumako'
