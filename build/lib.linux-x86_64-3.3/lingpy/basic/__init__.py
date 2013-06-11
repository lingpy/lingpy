# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-08 09:33
# modified : 2013-04-04 22:48
"""
This module provides basic classes for the handling of linguistic data.

The basic idea is to provide classes that allow the user to handle basic
linguistic datatypes (spreadsheets, wordlists) in a consistent way.

"""

__author__="Johann-Mattis List"
__date__="2013-04-04"

# set the general path to lingpy 
import os
_abs_path = os.path.split(
        os.path.abspath(
            __file__
            )
        )[0].replace('basic','')

# add imorts for basic classes
from .wordlist import Wordlist
from .spreadsheet import Spreadsheet
from .qlcwordlist import QLCWordlist
