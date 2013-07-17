# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-17 10:40
# modified : 2013-07-17 10:40
"""
Module handels all global parameters used in a LingPy session.
"""

__author__="Johann-Mattis List"
__date__="2013-07-17"

# builtin imports
from datetime import datetime,date
import os

# set the absolute path to lingpy
_abs_path = os.path.split(
        os.path.abspath(
            __file__
            )
        )[0]

# dictionary stores basic parameters that are used during a LingPy session
rcParams = dict(
        verbose = False,
        file_written               = "[i] Data has been written to file <{0}>.",
        missing_module             = "[WARNING] Module {0} could not be loaded. Some methods may not work properly.",
        filename                   = 'lingpy-'+str(date.today()),
        timestamp                  = str(datetime.today()),
        _path = _abs_path,
        input_error = "[!] Wrong input format, make sure the input format is of type {0}."
        )


