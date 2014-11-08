# *-* coding: utf-8 *-*
"""Implements the lingpy cache.

Some operations in lingpy may be time consuming, so we provide a mechanism to cache the
results of these operations.
"""
from __future__ import unicode_literals, print_function, absolute_import, division
import pickle

from pathlib import Path
from appdirs import user_cache_dir


DIR = Path(user_cache_dir('lingpy'))


def path(filename):
    if not isinstance(filename, Path):
        filename = Path(filename)
    return DIR.joinpath(filename.name + '.pkl')


def load(filename):
    with path(filename).open('rb') as fp:
        return pickle.load(fp)


def dump(data, filename):
    if not DIR.exists():
        DIR.mkdir(parents=True)  # pragma: no cover
    with path(filename).open('wb') as fp:
        pickle.dump(data, fp)
