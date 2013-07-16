# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-07-10 11:23
# modified : 2013-07-10 11:23
"""
Package provides classes and functions for testing and debugging
"""

__author__="Johann-Mattis List"
__date__="2013-07-10"

from datetime import date,datetime

# create a global timestamp variable
def _timestamp(time=''):
    """
    Function returns a timestamp.
    """
    if time in ['today','']:
        return str(date.today())
    elif time in ['now']:
        return str(datetime.today())
