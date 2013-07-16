# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-01-28 11:49
# modified : 2013-01-28 11:49
"""
Module provides common warnings for LingPy.
"""

__author__="Johann-Mattis List"
__date__="2013-01-28"

class ThirdPartyModuleError(Exception):
    """
    Common exception to be raised if third party modules are not installed.
    """
    def __init__(self,value):
        self.value = "Import of third party module <{0}> failed.".format(value)
        self.value += " Some methods won't work."

    def __str__(self):
        return repr(self.value)

    def warning(self):
        
        return "[WARNING] "+self.value

class InputFileError(Exception):
    """
    Basic exception to be raised if the input file contains obvious errors.
    """
    def __init__(self,value):
        self.value = "Erroneous coding of input file <{0}>. If it contains an ID column, make sure it consists only of integers.".format(value)

    def warning(self):

        return "[WARNING] "+self.value
