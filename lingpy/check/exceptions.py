# created: 2013-01-25
# modified: 2013-01-25
"""
Module provides common warnings for LingPy.
"""
__author__ = "Johann-Mattis List"
__date__ = "2013-01-25"

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
        
        print("[WARNING] "+self.value)
        return
