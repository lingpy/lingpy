# created :  2012-12-30
# modified : 2012-12-30

"This module provides basic functions and classes."

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
