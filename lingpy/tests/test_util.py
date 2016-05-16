from __future__ import unicode_literals
from unittest import TestCase

from lingpy.tests.util import WithTempDir, test_data
from lingpy import util


class Test(WithTempDir):
    def test_write_text_file(self):
        def lines_generator(n):
            for i in range(n):
                yield 'line%s' % i

        path = self.tmp_path('test')
        util.write_text_file(path, 'test')
        self.assertEqual(util.read_text_file(path), 'test')

        util.write_text_file(path, ['line1', 'line2'])
        self.assertEqual(len(util.read_text_file(path, lines=True)), 2)

        util.write_text_file(path, lines_generator(5))
        self.assertEqual(len(util.read_text_file(path, lines=True)), 5)

    def test_TextFile(self):
        path = self.tmp_path('test')
        with util.TextFile(path) as fp:
            fp.writelines(['line1\n', 'line2\n'])
        self.assertEqual(len(util.read_text_file(path, lines=True)), 2)


class TestCombinations(TestCase):
    def test_combinations2(self):
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

        for l in [list(range(5)), 'abcdefg']:
            self.assertEqual(list(util.combinations2(l)), list(f(l)))
            self.assertEqual(list(util.multicombinations2(l)), list(fm(l)))


class TestJoin(TestCase):
    def test_join(self):
        self.assertEqual(util.join('.'), '')
        self.assertEqual(util.join('.', 1), '1')
        self.assertEqual(util.join('.', 1, 2), '1.2')

    def test_dotjoin(self):
        self.assertEqual(util.dotjoin(1, 2), '1.2')
        self.assertEqual(util.dotjoin([1, 2]), '1.2')
        self.assertEqual(util.dotjoin((1, 2)), '1.2')
        self.assertEqual(
            util.dotjoin((i for i in range(1, 3)), condition=lambda j: j > 1), '2')
        self.assertEqual(util.dotjoin(i for i in range(1, 3)), '1.2')

def test_as_string():

    out = util.as_string('text', pprint=False)
    assert out == 'text' 

def test_read_csv_file():

    lines1 = util.read_csv_file(test_data('mycsvwordlist.csv'))
