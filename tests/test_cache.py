from lingpy import cache


def test_cache(tmppath):
    d = {'a': 123}
    filename = 'lingpy_test.CSV'
    cache.dump(d, filename, d=tmppath / 'cache')
    assert cache.load(filename, d=tmppath / 'cache') == d
