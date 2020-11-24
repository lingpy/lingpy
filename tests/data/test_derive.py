import pytest

from lingpy import cache
from lingpy.data.derive import compile_model, compile_dvt


@pytest.fixture
def cache_dir(tmppath):
    d = tmppath / '_test'
    d.mkdir()
    d.joinpath('converter').write_text("""\
p : p, ɸ, p͡f
b : b, β, b͡v
f : f
v : v
m : m, ɱ
w : w, ɰ, ʋ, ʍ
8 : θ, θ, Ɵ, ð""", encoding='utf8')
    d.joinpath('scorer').write_text("""\
p : c, b:1, f:2
b : c, -
f : c, -
v : c, -
m : v, w:1
w : v, m:1
8 : t, -""")
    return d.parent


def test_compile_model(cache_dir):
    compile_model('_test', str(cache_dir))
    sound_classes = cache.load('_test.converter')
    assert sound_classes['b'] == 'b'
    assert cache_dir.joinpath('_test', 'matrix').exists()


def test_compile_dvt():
    compile_dvt()
    assert len(cache.load('dvt')) == 3
