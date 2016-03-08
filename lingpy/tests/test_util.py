from __future__ import unicode_literals

from lingpy.tests.util import WithTempDir
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
