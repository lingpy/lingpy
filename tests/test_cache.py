from lingpy import cache


def test_cache(tmp_path):
    d = {'a': 123}
    filename = 'lingpy_test.CSV'
    cache.dump(d, filename, d=tmp_path / 'cache')
    assert cache.load(filename, d=tmp_path / 'cache') == d
