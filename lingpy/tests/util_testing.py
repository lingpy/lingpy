from __future__ import unicode_literals

import sys
import tempfile
import unittest
import contextlib

from six import StringIO

from clldutils.path import Path, rmtree


class WithTempDirMixin(object):
    """Composable test fixture providing access to a temporary directory.

    http://nedbatchelder.com/blog/201210/multiple_inheritance_is_hard.html
    """

    def setUp(self):
        super(WithTempDirMixin, self).setUp()
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        rmtree(self.tmp, ignore_errors=True)
        super(WithTempDirMixin, self).tearDown()

    def tmp_path(self, *comps):
        return self.tmp.joinpath(*comps)


class WithTempDir(WithTempDirMixin, unittest.TestCase):
    """Backwards compatible test base class."""


@contextlib.contextmanager
def capture(func, *args, **kw):
    with capture_all(func, *args, **kw) as res:
        yield res[1]


@contextlib.contextmanager
def capture_all(func, *args, **kw):
    out, sys.stdout = sys.stdout, StringIO()
    err, sys.stderr = sys.stderr, StringIO()
    ret = func(*args, **kw)
    sys.stdout.seek(0)
    sys.stderr.seek(0)
    yield ret, sys.stdout.read(), sys.stderr.read()
    sys.stdout, sys.stderr = out, err
