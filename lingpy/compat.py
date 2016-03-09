"""Functionality to provide compatibility across supported python versions"""
from __future__ import unicode_literals, division, absolute_import, print_function

# python 2.7 compatibility: see http://stackoverflow.com/a/21368622
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError
