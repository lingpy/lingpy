from lingpy import util


def test_write_text_file(tmp_path):
    def lines_generator(n):
        for i in range(n):
            yield 'line%s' % i

    path = tmp_path / 'test'
    util.write_text_file(path, 'test')
    assert util.read_text_file(path) == 'test'

    util.write_text_file(path, ['line1', 'line2'])
    assert len(util.read_text_file(path, lines=True)) == 2

    util.write_text_file(path, lines_generator(5))
    assert len(util.read_text_file(path, lines=True)) == 5


def test_TextFile(tmp_path):
    path = tmp_path / 'test'
    with util.TextFile(path) as fp:
        fp.writelines(['line1\n', 'line2\n'])
    assert len(util.read_text_file(path, lines=True)) == 2


def test_combinations2():
    def f(l):
        for i, a1 in enumerate(l):
            for j, a2 in enumerate(l):
                if i < j:
                    yield a1, a2

    def fm(l):
        for i, a1 in enumerate(l):
            for j, a2 in enumerate(l):
                if i <= j:
                    yield a1, a2

    for ch in [list(range(5)), 'abcdefg']:
        assert list(util.combinations2(ch)) == list(f(ch))
        assert list(util.multicombinations2(ch)) == list(fm(ch))


def test_join():
    assert util.join('.') == ''
    assert util.join('.', 1) == '1'
    assert util.join('.', 1, 2) == '1.2'


def test_dotjoin():
    assert util.dotjoin(1, 2) == '1.2'
    assert util.dotjoin([1, 2]) == '1.2'
    assert util.dotjoin((1, 2)) == '1.2'
    assert util.dotjoin((i for i in range(1, 3)), condition=lambda j: j > 1) == '2'
    assert util.dotjoin(i for i in range(1, 3)) == '1.2'


def test_as_string():
    out = util.as_string('text', pprint=False)
    assert out == 'text'
