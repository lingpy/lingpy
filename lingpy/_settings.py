# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-17 10:40
# modified : 2013-07-17 10:40
"""
Module handels all global parameters used in a LingPy session.
"""

__author__ = "Johann-Mattis List"
__date__ = "2013-07-17"

# builtin imports
from datetime import datetime, date
import os

# initialize rcParams with filename and timestamp
rcParams = dict(
    # set the absolute path to lingpy
    _path=os.path.split(os.path.abspath(__file__))[0],
    filename='lingpy-' + str(date.today()),
    timestamp=datetime.strftime(datetime.today(),'%Y-%m-%d %H:%M'),
    answer_yes=['y', 'Y', 'j', 'J', 'yes'],
    verbose=False,
    debug=False,
    basic_orthography='fuzzy',
)
