from lingpy.read.starling import star2qlc


def test_star2qlc(test_data):
    star2qlc(str(test_data / 'rom.starling.tsv'), debug=True)
    res = star2qlc(str(test_data / 'rom.starling.tsv'))
    assert len(res) == 208
