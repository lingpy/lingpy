"""
Module handels all global parameters used in a LingPy session.
"""
from __future__ import unicode_literals, division, print_function
from datetime import datetime, date


# initialize rcParams with filename and timestamp
rcParams = dict(
    filename='lingpy-' + str(date.today()),
    timestamp=datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M'),
    basic_orthography='fuzzy',
)
