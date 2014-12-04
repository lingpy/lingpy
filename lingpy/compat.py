"""Functionality to provide compatibility across supported python versions"""
from __future__ import unicode_literals, division, absolute_import, print_function
import json

# python 2.7 compatibility: see http://stackoverflow.com/a/21368622
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

from six import PY3


def jsonload(path):
    """python 2 + 3 compatible version of json.load.

    :return: The python object read from path.
    """
    kw = {}
    if PY3:  # pragma: no cover
        kw['encoding'] = 'utf8'
    with open(path, **kw) as fp:
        return json.load(fp)
