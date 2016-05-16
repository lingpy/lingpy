# *-* coding: utf-8 *-*
from __future__ import unicode_literals, print_function, division, absolute_import

from lingpy.tests.util import WithTempDir
from lingpy import cache


SCORER = """\
p : c, b:1, f:2
b : c, -
f : c, -
v : c, -
m : v, w:1
w : v, m:1
8 : t, -"""


class TestDerive(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.tmp_path('_test').mkdir()
        with self.tmp_path('_test', 'converter').open(mode='w', encoding='utf8') as fp:
            fp.write("""\
p : p, ɸ, p͡f
b : b, β, b͡v
f : f
v : v
m : m, ɱ
w : w, ɰ, ʋ, ʍ
8 : θ, θ, Ɵ, ð""")

        with self.tmp_path('_test', 'scorer').open(mode='w', encoding='utf8') as fp:
            fp.write(SCORER)

    def test_compile_model(self):
        from lingpy.data.derive import compile_model

        compile_model('_test', self.tmp)
        sound_classes = cache.load('_test.converter')
        self.assertEqual(sound_classes['b'], 'b')
        self.assertTrue(self.tmp_path('_test', 'matrix').exists())

    def test_compile_dvt(self):
        from lingpy.data.derive import compile_dvt

        compile_dvt()
        self.assertEqual(len(cache.load('dvt')), 3)
