"""Utilities used in lingpy tests"""
import os


def test_data(*comps):
    """Access test data files.

    :param comps: Path components of the data file path relative to the test_data dir.
    :return: Absolute path to the specified test data file.
    """
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'test_data', *comps)
