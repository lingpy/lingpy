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
