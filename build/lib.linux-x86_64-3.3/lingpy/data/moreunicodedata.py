"""
This module extend the Python standard library unicodedata module by providing access to Unicode scripts.txt and blocks.txt.
"""

__author__ = "Steven Moran"
__date__="2013-01"

from unicodedata import *
import urllib.request
import regex as re

class MoreUnicodedata:
    """
    A class that extends the functionality of the standard Python libary's unicodedata module that includes
    accessors to scripts and blocks information.
    """

    def __init__(self):
        self.scripts_url = "http://www.unicode.org/Public/UNIDATA/Scripts.txt"
        self.blocks_url = "http://www.unicode.org/Public/UNIDATA/Blocks.txt"
        self.script_data = {}
        self.block_data = {}
        self._compile_scripts_data()
        self._compile_blocks_data()

    def _compile_blocks_data(self):
        """
        Process Unicode Blocks.txt and store the data.
        """
        ids = []
        blocks = []
        file = urllib.request.urlopen(self.blocks_url)

        for line in file:
            line = line.decode("utf-8") # urllib.request returns byte object; convert it to str
            p = re.findall('([0-9A-F]+)(?:\.\.([0-9A-F]+))\;\s(\w+.+)', line)
            if p:
                begin, end, block = p[0]

                if block not in blocks:
                    blocks.append(block)

                ids.append((int(begin, 16), int(end, 16), blocks.index(block)))
        ids.sort()

        self.block_data["blocks"] = blocks
        self.block_data["ids"] = ids

    def get_block(self, chr):
        """
        Given a Unicode character, returns the block it is assigned to.
        """
        l = 0
        r = len(self.block_data['ids']) - 1
        c = ord(chr)

        # if out of range...

        while r >= l:
            m = (l + r) >> 1
            if c < self.block_data['ids'][m][0]:
                r = m - 1
            elif c > self.block_data['ids'][m][1]:
                l = m + 1
            else:
                return self.block_data['blocks'][self.block_data['ids'][m][2]]
        return 'Unknown', 'Zzzz'        


    def _compile_scripts_data(self):
        """
        Process Unicode Scripts.txt and store the data.
        """
        ids = []
        scripts = []
        categories = []        

        file = urllib.request.urlopen(self.scripts_url) # get file

        for line in file:
            line = line.decode("utf-8") # urllib.request returns byte object; convert it to str
            p = re.findall('([0-9A-F]+)(?:\.\.([0-9A-F]+))?\W+(\w+)\s*#\s*(\w+)', line)
            if p:
                begin, end, script, category = p[0]
                
                if script not in scripts:
                    scripts.append(script)

                if category not in categories:
                    categories.append(category)
                    
                ids.append((int(begin, 16), int(end or begin, 16), scripts.index(script), categories.index(category)))
        ids.sort()

        self.script_data["scripts"] = scripts
        self.script_data["categories"] = categories
        self.script_data["ids"] = ids


    def get_script_category(self, chr):
        """
        Given a Unicode character, returns a tuple of (scriptname, category).
        """
        l = 0
        r = len(self.script_data['ids']) - 1
        c = ord(chr)

        # if out of range...

        while r >= l:
            m = (l + r) >> 1
            if c < self.script_data['ids'][m][0]:
                r = m - 1
            elif c > self.script_data['ids'][m][1]:
                l = m + 1
            else:
                return (
                    self.script_data['scripts'][self.script_data['ids'][m][2]], 
                    self.script_data['categories'][self.script_data['ids'][m][3]])
        return 'Unknown', 'Zzzz'        

    def get_script(self, chr):
        """
        Returns the script for a given Unicode character.
        """
        script, _ = self.get_script_category(chr)
        return script

    def get_category(self, chr):
        """
        Returns the category for a given Unicode character.
        """
        _, category = self.get_script_category(chr)
        return category


if __name__=="__main__":
    u = MoreUnicodedata()
    print(u.get_block("a"))
    # print(u.get_script("a"))
    
