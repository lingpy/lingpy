from lingpy.basic.wordlist import from_cldf
from lingpy.convert.cldf import to_cldf


def test_from_cldf(test_data, tmp_path):
    wl = from_cldf(str(test_data / 'cldf/test-metadata.json'), language='Name',
                   concept='Name', concepticon="Concepticon_ID",
                   glottocode='glottocode')

    assert wl.width == 29
    assert wl.height == 1
    assert wl.entries[0] == 'alignment'
    assert wl.cols[0] == 'Anuta'
    assert wl.cols[28] == 'Vaeakau-Taumako'

    to_cldf(wl, path=tmp_path)
    assert tmp_path.joinpath('Wordlist-metadata.json').exists()
