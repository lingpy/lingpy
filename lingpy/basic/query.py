# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-01-27 16:46
# modified : 2014-01-27 16:46
"""
A wrapper module which can be used to query wordlists.
"""

__author__="Johann-Mattis List"
__date__="2014-01-27"

import os

from .wordlist import Wordlist
from ..settings import rc, rcParams


class Query(Wordlist):

    def __init__(self, filename, row="concept", col="doculect", conf="",
            settings=None):

        Wordlist.__init__(filename, row, col, conf)

    def select(



