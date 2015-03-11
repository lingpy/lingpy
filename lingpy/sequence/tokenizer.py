# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
"""
This module provides graphemic and orthographic tokenization using with orthography profiles.
"""

__author__ = "Steven Moran"
__date__ = "2010-12-01"

import os
import logging
import codecs
import unicodedata

# basic lingpy imports
from ..settings import rcParams
from .. import log
from .. import util

try:
    import regex as re
except ImportError:
    import re
    log.missing_module('regex')


class Tokenizer(object):
    """
    Class for orthographic parsing using orthography profiles as designed for the QLC project.

    Parameters
    ----------

    orthography_profile : string (default = None)
        Filename of the a document source-specific orthography profile and rules file.

    Notes
    -----

    The Tokenizer reads in an orthography profile and calls a helper 
    class to build a tree data structure, which stores the possible Unicode 
    character combinations that are specified in the orthography profile 
    and appear in the data source.

    For example, an orthography profile might specify that in source X 
    <uu> is a single grapheme (Unicode parlance: tailored grapheme) and 
    therefore it should be chunked as so. Given an orthography profile and 
    some data to parse, the process would look like this:

    input string example: uubo uubo
    output string example: uu b o # uu b o

    See also the tokenizer examples in lingpy/scripts/tokenize

    Additionally, the Tokenizer provides functionality to transform graphemes 
    into associated character(s) specified in additional columns in the orthography 
    profile. A dictionary is created that keeps a mapping between source-specific 
    graphemes and their counterparts (e.g. an IPA column in the orthography profile).

    The tokenizer can also be used for pure Unicode character and grapheme 
    tokenization, i.e. it uses the Unicode standard grapheme parsing rules, as 
    implemented in the Python regex package by Matthew Barnett, to do basic tokenization 
    with the "\X" grapheme regular expression match. This grapheme match 
    combines one or more Combining Diacritical Marks to their base character. 
    These are called "Grapheme clusters" in Unicode parlance. With these functions 
    the Tokenizer is meant to do basic rudimentary parsing for things like getting 
    an additional unigram model (segments and their counts) in an input data source.

    An additional method (in its infancy) called combine_modifiers handles the 
    case where there are Unicode Spacing Modifier Letters, which are not explicitly 
    combined to their base character in the Unicode Standard. These graphemes 
    are called "Tailored grapheme clusters" in Unicode. For more information 
    see the Unicode Standard Annex #29: Unicode Text Segmentation:

    http://www.unicode.org/reports/tr29/

    Lastly, the Tokenizer can be used to transformation as specified in an 
    orthography rules file. These transformations are specified in a separate 
    file from the orthography profile (that specifics the document specific graphemes, 
    and possibly their IPA counterparts) and the orthography rules should 
    be applied to the output of an OrthographyParser.

    In an orthography rules file, rules are given in order in regular 
    expressions, e.g. this rule replaces a vowel followed by an <n> 
    followed by <space> followed by a second vowel with first vowel 
    <space> <n> <space> second vowel, e.g.:

    ([a|á|e|é|i|í|o|ó|u|ú])(n)(\s)([a|á|e|é|i|í|o|ó|u|ú]), \1 \2 \4

    Examples
    --------

    >>> from lingpy.sequence.tokenizer import *
    >>> t = Tokenizer("test.prf")
    >>> word = "baach"
    >>> print(t.characters(word))
    b a a c h
    >>> print(t.graphemes(word))
    b aa ch
	>>> print(t.transform(word, "ipa"))
	b aː tʃ    

    """
    def __init__(self, orthography_profile=None):
        if orthography_profile and not os.path.exists(orthography_profile):
            raise ValueError("The orthography profile you specified does not exist!")
        self.orthography_profile = orthography_profile
        self.orthography_profile_rules = None
        self.column_labels = None

        # orthography profile processing
        if self.orthography_profile:
            # read in orthography profile and create a tree structure for tokenization
            self.root = createTree(self.orthography_profile)

            # store column labels from the orthography profile
            self.column_labels = []

            # look up table of graphemes to other column transforms
            self.mappings = {}

            # double check that there are no duplicate graphemes in the orthography profile
            self.op_graphemes = {}

            # process the orthography profiles and rules
            self._init_profile(self.orthography_profile)

            rules_path = os.path.splitext(self.orthography_profile)[0] + '.rules'
            if os.path.isfile(rules_path):
                self.orthography_profile_rules = rules_path
                self.op_rules = []
                self.op_replacements = []
                self._init_rules(self.orthography_profile_rules)
        else:
            try:
                import regex as re
            except ImportError:
                raise ImportError(
                    "Please install the `regex` module to use Tokenizer without an orthography_profile."
                )
            self.grapheme_pattern = re.compile("\X", re.UNICODE)
            

        log.debug("Orthography profile: %s" % self.orthography_profile)
        log.debug("Orthography rules: %s" % self.orthography_profile_rules)
        log.debug("Columns labels: %s" % self.column_labels)

    def _init_profile(self, f):
        # Process and initialize data structures given an orthography profile.
        for line_count, line in enumerate(util.read_config_file(f, normalize='NFD')):
            # deal with the columns header -- should always start with "graphemes" as per
            # the orthography profiles specification
            if line.lower().startswith("graphemes"):
                column_tokens = line.split("\t")

                # clean the header
                for column_token in column_tokens:
                    self.column_labels.append(column_token.lower().strip())
                continue

            # split the orthography profile into columns
            tokens = line.split("\t") 
            grapheme = tokens[0].strip()

            # check for duplicates in the orthography profile (fail if dups)
            if not grapheme in self.op_graphemes:
                self.op_graphemes[grapheme] = 1
            else:
                raise Exception("You have a duplicate in your orthography profile.")

            if len(tokens) == 1:
                continue

            for i, token in enumerate(tokens):
                token = token.strip()
                self.mappings[grapheme, self.column_labels[i].lower()] = token
                log.debug('%s %s' % (grapheme, self.column_labels[i].lower()))

        # print the tree structure if debug mode is on
        if log.get_logger().getEffectiveLevel() <= logging.INFO:
            log.debug("A graphical representation of your orthography profile in a tree ('*' denotes sentinels):\n")
            printTree(self.root, "")
            print()

    def _init_rules(self, f):
        # Process the orthography rules file.
        for line in util.read_config_file(f, normalize='NFD'):
            rule, replacement = line.split(",")
            rule = rule.strip()  # just in case there's trailing whitespace
            replacement = replacement.strip()  # because there's probably trailing whitespace!
            self.op_rules.append(re.compile(rule))
            self.op_replacements.append(replacement)

        # check that num rules == num replacements; if not fail
        if len(self.op_rules) != len(self.op_replacements):
            raise ValueError("Number of inputs does not match number of outputs in the rules file.")

    def characters(self, string):
        """
        Given a string as input, return a space-delimited string of Unicode characters (code points rendered as glyphs).

        Parameters
        ----------
        string : str
            A Unicode string to be parsed into graphemes.

        Returns
        -------
        result : str
            String returned is space-delimited on Unicode characters and contains "#" to mark word boundaries.
            The string is in NFD.

        Notes
        -----
        Input is first normalized according to Normalization Ford D(ecomposition).
        String returned contains "#" to mark word boundaries.
        """

        string = string.replace(" ", "#") # add boundaries between words
        string = unicodedata.normalize("NFD", string)
        result = ""
        for character in string:
            result += character+" "
        return result.strip()

    def grapheme_clusters(self, string):
        """
        Given a string as input, return a space-delimited string of Unicode graphemes using the "\X" regular expression.

        Parameters
        ----------
        string : str
            A Unicode string to be parsed into graphemes.

        Returns
        -------
        result : str
            String returned is space-delimited on Unicode graphemes and contains "#" to mark word boundaries.
            The string is in NFD.

        Notes
        -----
        Input is first normalized according to Normalization Ford D(ecomposition).

        See: Unicode Standard Annex #29: UNICODE TEXT SEGMENTATION
        http://www.unicode.org/reports/tr29/
        """

        # init the regex Unicode grapheme cluster match
        return ' '.join(self.grapheme_pattern.findall(
            unicodedata.normalize("NFD", string.replace(" ", "#"))))


    def graphemes(self, string):
        """
        Tokenizes strings given an orthograhy profile that specifies graphemes in a source doculect.

        Parameters
        ----------
        string : str
            The str to be parsed and formatted.

        Returns
        -------
        result : str
            The result of the parsed and QLC formatted str.

        """
        string = unicodedata.normalize("NFD", string)
        
        # if no orthography profile is specified, simply return 
        # Unicode grapheme clusters, regex pattern "\X"
        if self.orthography_profile == None:
            return self.grapheme_clusters(string)

        parses = []
        for word in string.split():
            parse = getParse(self.root, word)

            # case where the parsing fails
            if len(parse) == 0:
                # replace characters in string but not in orthography profile with <?>
                parse = " "+self.find_missing_characters(self.characters(word))
                # write problematic stuff to standard error
                log.debug("The string '{0}' does not parse given the specified orthography profile {1}.\n".format(word, self.orthography_profile))
            
            parses.append(parse)

        # remove the outter word boundaries
        result = "".join(parses).replace("##", "#")
        result = result.rstrip("#")
        result = result.lstrip("#")
        return result.strip()


    def transform(self, string, column="graphemes"):
        """
        Transform a string's graphemes into the mappings given in a different column 
        in the orthography profile. By default this function returns an orthography 
        profile grapheme tokenized string.

        Parameters
        ----------
        string : str
            The input string to be parsed.

        conversion : str (default = "graphemes")
            The label of the column to transform to. Default it to tokenize with orthography profile.

        Returns
        -------
        result : str
            Result of the transformation.

        """
        # column labels are normalized 
        column = column.lower()

        # This method can't be called unless an orthography profile was specified.
        if not self.orthography_profile:
            raise Exception("This method only works when an orthography profile is specified.")

        if column == "graphemes":
            return self.graphemes(string)

        # if the column label for conversion doesn't exist, return grapheme tokenization
        if column not in self.column_labels:
            return self.graphemes(string)

        # first tokenize the input string into orthography profile graphemes
        tokenized_string = self.graphemes(string)
        tokens = tokenized_string.split()

        result = []
        for token in tokens:
            # special cases: word breaks and unparsables
            if token == "#":
                result.append("#")
                continue
            if token == "?":
                result.append("?")
                continue

            # transform given the grapheme and column label; skip NULL
            target = self.mappings[token, column]
            if not target == "NULL":
                # result.append(self.mappings[token, column]) 
                result.append(target)

        return " ".join(result).strip()

    def tokenize(self, string, column="graphemes"):
        """
        This function determines what to do given any combination 
        of orthography profiles and rules or not orthography profiles
        or rules.

        Parameters
        ----------
        string : str
            The input string to be tokenized.

        column : str (default = "graphemes")
            The column label for the transformation, if specified.

        Returns
        -------
        result : str
            Result of the tokenization.

        """

        # column labels are normalized 
        column = column.lower() 

        if self.orthography_profile and self.orthography_profile_rules:
            return self.rules(self.transform(string, column))

        if not self.orthography_profile and not self.orthography_profile_rules:
            return self.grapheme_clusters(string)

        if self.orthography_profile and not self.orthography_profile_rules:
            return self.transform(string, column)

        # it's not yet clear what the order for this procedure should be
        if not self.orthography_profile and self.orthography_profile_rules:
            return self.rules(self.grapheme_clusters(string))

    def transform_rules(self, string):
        """
        Convenience function that first tokenizes a string into orthographic profile-
        specified graphemes and then applies the orthography profile rules.

        Parameters
        ----------
        string : str
            The input string to be transformed.

        Returns
        -------
        result : str
            Result of the transformation.

        """
        return self.rules(self.transform(string))

    def rules(self, string):
        """
        Function to parse input string and return output of str with ortho rules applied.

        Parameters
        ----------
        string : str
            The input string to be parsed.

        Returns
        -------
        result : str
            Result of the orthography rules applied to the input str.

        """
        # if no orthography profile was initiated, this method can't be called
        # if self.orthography_profile == None:
        #    raise Exception("This function requires that an orthography profile is specified.")

        # if no orthography profile rules file has been specified, simply return the string
        if self.orthography_profile_rules == None:
            return string

        result = unicodedata.normalize("NFD", string)
        for i in range(0, len(self.op_rules)):
            match = self.op_rules[i].search(result)
            if match:
                result = re.sub(self.op_rules[i], self.op_replacements[i], result)
                log.debug("Input/output:"+"\t"+string+"\t"+result)
                log.debug("Pattern/replacement:"+"\t"+self.op_rules[i].pattern+"\t"+self.op_replacements[i])

        # this is incase someone introduces a non-NFD ordered sequence of characters
        # in the orthography profile
        result = unicodedata.normalize("NFD", result)
        return result

    def find_missing_characters(self, char_tokenized_string):
        """
        Given a string tokenized into characters, return a characters
        tokenized string where each character missing from the orthography 
        profile is replaced with a question mark <?>.

        Parameters
        ----------
        string : str
            A character tokenized string.

        Returns
        -------
        result : str
            Result of the tokenization.

        """
        result = []
        chars = char_tokenized_string.split()
        for char in chars:
            if not char in self.op_graphemes:
                result.append("?")
            else:
                result.append(char)
        return " ".join(result)

    def tokenize_ipa(self, string):
        # Experimental method for tokenizing IPA.
        return self.combine_modifiers(self.grapheme_clusters(string))

    def combine_modifiers(self, string):
        """
        Given a string that is space-delimited on Unicode grapheme cluters, 
        group Unicode modifier letters with their preceeding base characters.

        Parameters
        ----------
        string : str
            A Unicode string tokenized into grapheme clusters to be tokenized into simple IPA.

        .. todo:: check if we need to apply NDF after string is parsed

        """
        result = []
        graphemes = string.split()
        temp = ""
        count = len(graphemes)

        for grapheme in reversed(graphemes):
            count -= 1
            if len(grapheme) == 1 and unicodedata.category(grapheme) == "Lm":
                temp = grapheme + temp
                # hack for the cases where a space modifier is the first character in the str
                if count == 0:
                    result[-1] = temp + result[-1]
                continue

            result.append(grapheme + temp)
            temp = ""

        # check for tie bars
        segments = result[::-1]

        i = 0
        r = []
        while i < len(segments):
            if ord(segments[i][-1]) in [865, 860]:
                r.append(segments[i]+segments[i+1])
                i = i+2
            else:
                r.append(segments[i])
                i += 1

        return " ".join(r)

    def exists_multiple_columns(self):
        """
        Returns boolean of whether multiple columns exist in the orthography profile, e.g. graphemes and IPA and x, etc.
        """
        if len(self.column_labels) > 1:
            return True
        else:
            return False

    def remove_spaces(self, string):
        string = string.lstrip("# ")
        string = string.rstrip(" #")
        string = re.sub("\s+", "", string)
        string = string.replace("#", " ")
        return string


# ---------- Tree node --------
class TreeNode(object):
    """
    Private class that creates the tree data structure from the orthography profile for parsing.
    """

    def __init__(self, char):
        self.char = char
        self.children = {}
        self.sentinel = False

    def isSentinel(self):
        return self.sentinel

    def getChar(self):
        return self.char

    def makeSentinel(self):
        self.sentinel = True

    def addChild(self, char):
        child = self.getChild(char)
        if not child:
            child = TreeNode(char)
            self.children[char] = child
        return child

    def getChild(self, char):
        if char in self.children:
            return self.children[char]
        else:
            return None

    def getChildren(self):
        return self.children


# ---------- Util functions ------
def createTree(file_name):
    # Internal function to add a multigraph starting at node.
    def addMultigraph(node, line):
        for char in line:
            node = node.addChild(char)
        node.makeSentinel()

    # Add all multigraphs in each line of file_name. Skip "#" comments and blank lines.
    root = TreeNode('')
    root.makeSentinel()

    f = codecs.open(file_name, "r", "utf-8")
    header = []

    for line in f:
        line = line.strip()

        # skip any comments
        if line.startswith("#") or line == "":
            continue

        # deal with the columns header -- should always start with "graphemes" as per the orthography profiles specification
        if line.lower().startswith("graphemes"):
            header = line.split("\t")
            continue

        line = unicodedata.normalize("NFD", line)
        tokens = line.split("\t") # split the orthography profile into columns
        grapheme = tokens[0]
        addMultigraph(root, grapheme)
    f.close()
    return root


def printMultigraphs(root, line, result):
    # Base (or degenerate..) case.
    if len(line) == 0:
        result += "#"
        return result

    # Walk until we run out of either nodes or characters.
    curr = 0   # Current index in line.
    last = 0   # Index of last character of last-seen multigraph.
    node = root
    while curr < len(line):
        node = node.getChild(line[curr])
        if not node:
            break
        if node.isSentinel():
            last = curr
        curr += 1

    # Print everything up to the last-seen sentinel, and process
    # the rest of the line, while there is any remaining.
    last = last + 1  # End of span (noninclusive).
    result += line[:last]+" "
    return printMultigraphs(root, line[last:], result)


def getParse(root, line):
    parse = getParseInternal(root, line)
    if len(parse) == 0:
        return ""
    return "# " + parse


def getParseInternal(root, line):
    # Base (or degenerate..) case.
    if len(line) == 0:
        return "#"

    parse = ""
    curr = 0
    node = root
    while curr < len(line):
        node = node.getChild(line[curr])
        curr += 1
        if not node:
            break
        if node.isSentinel():
            subparse = getParseInternal(root, line[curr:])
            if len(subparse) > 0:
                # Always keep the latest valid parse, which will be
                # the longest-matched (greedy match) graphemes.
                parse = line[:curr] + " " + subparse

    # Note that if we've reached EOL, but not end of valid grapheme,
    # this will be an empty string.
    return parse


def printTree(root, path):
    for char, child in root.getChildren().items():
        if child.isSentinel():
            char += "*"
        branch = (" -- " if len(path) > 0 else "")
        printTree(child, path + branch + char)
    if len(root.getChildren()) == 0:
        print(path)
