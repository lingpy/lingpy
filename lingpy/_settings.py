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

# initialize rcParams with filename and timestamp
rcParams = dict(
    filename='lingpy-' + str(date.today()),
    timestamp=datetime.strftime(datetime.today(),'%Y-%m-%d %H:%M'),
    answer_yes=['y', 'Y', 'j', 'J', 'yes'],
    basic_orthography='fuzzy',
)
