# created: 2013-01-25
# modified: 2013-01-25
"""
Module provides basic messages that are frequently used throughout LingPy.
"""
__author__ = "Johann-Mattis List"
__date__ = "2013-01-25"

class FileWriteMessage(object):

    def __init__(self,filename,fileformat):

        self.value = filename+'.'+fileformat
    
    def message(self,mtype):

        if mtype == 'written':
            out = "[i] Data has been written to file <{0}>.".format(
                    self.value
                    )

        print(out)
        
        return

    def __str__(self):
        return "[i] Data has been written to file <{0}>.".format(self.value)

class LoadDataMessage(object):

    def __init__(self,*args):
        
        if len(args) > 1:
            self.value = ', '.join(args[:-1])+', and '+args[-1]
        else:
            self.value = args[0]

    def message(self,mtype):
        
        if mtype == 'loaded':
            out = "[i] Successfully loaded {0}.".format(self.value)
        elif mtype == 'failed':
            out = "[i] Failed to load {0}.".format(self.value)

        print(out)
        return

class LingPyDeprecationWarning(DeprecationWarning):
    """
    Class serves as a basis for deprecation warnings.
    """

    def __init__(self,*args):

        self.args = args

    def __str__(self):

        if len(self.args) == 2:
            return "[WARNING] Use of '{0[0]}' is deprecated. Use '{0[1]}' instead.".format(self.args)
        else:
            return str(self.args)

