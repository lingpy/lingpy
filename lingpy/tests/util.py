"""Utilities used in lingpy tests"""
import os
import unittest
from tempfile import mkdtemp
import shutil

from pathlib import Path

from lingpy.util import lingpy_path


def test_data(*comps):
    """Access test data files.

    :param comps: Path components of the data file path relative to the test_data dir.
    :return: Absolute path to the specified test data file.
    """
    return lingpy_path('tests', 'test_data', *comps)


class WithTempDir(unittest.TestCase):
    def setUp(self):
        self.tmp = mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def tmp_path(self, *comps):
        return Path(self.tmp).joinpath(*comps)
