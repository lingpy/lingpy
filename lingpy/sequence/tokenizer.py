"""
This module provides graphemic and orthographic tokenization using with orthography profiles.
"""

__author__ = "Steven Moran"
__date__ = "2010-12-01"

import sys
import codecs
import unicodedata

# basic lingpy imports
from ..settings import rcParams

try:
    import regex as re
except ImportError:
    import re
    print(rcParams['W_missing_module'].format('regex'))

class Tokenizer(object):
    """
    Class for orthographic parsing using orthography profiles as designed for the QLC project.

    Parameters
    ----------

    orthography_profile : sting
        Filepath to a document source-specific orthography profile (e.g. test.prf).

    orthography_profile_rules : string
        Filepath to the orthography profile rules file (e.g. test.rules).

    debug : int (default = 0)
        Turn debug mode on / off.


    Notes
    -----

    The Tokenizer reads in an orthography profile and calls a helper 
    class to build a trie data structure, which stores the possible Unicode 
    character combinations that are specified in the orthography profile 
    and appear in the data source.

    For example, an orthography profile might specify that in source X 
    <uu> is a single grapheme (Unicode parlance: tailored grapheme) and 
    thererfore it should be chunked as so. Given an orthography profile and 
    some data to parse, the process would look like this:

    input string example: uubo uubo
    output string example: uu b o # uu b o

    See also the test_orthography script in lingpy/scripts/orthography/

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

    """


    def __init__(self, orthography_profile=None, orthography_profile_rules=None, debug=0):
        # path and filename of orthography profile
        self.orthography_profile = orthography_profile
        self.orthography_profile_rules = orthography_profile_rules

        self.debug = debug

        if self.orthography_profile:
            # read in orthography profile and create a trie structure for tokenization
            self.root = createTree(self.orthography_profile)

            # store column labels from the orthography profile
            self.column_labels = []

            # look up table of graphemes to other column transforms
            self.mappings = {}

            # double check that there are no duplicate graphemes in the orthography profile
            self.op_graphemes = {}

            # process the orthography profiles and rules
            self._init_profile(self.orthography_profile)

        # orthography profile rules and replacements
        if self.orthography_profile_rules:
            self.op_rules = []
            self.op_replacements = []
            self._init_rules(self.orthography_profile_rules)

    
    def _init_profile(self, file):
        """
        Process and initialize data structures given an orthography profile.
        """
        try:
            codecs.open(file, 'r', 'utf-8')
        except IOError as e:
            print("\n[i] There is no file at the path you've specified.\n\n")

        file = codecs.open(file, "r", "utf-8")
        line_count = 0
        for line in file:
            line_count += 1
            line = line.strip()

            # skip any comments in the orthography profile
            if line.startswith("#") or line == "":
                continue

            # deal with the columns header -- should always start with "graphemes" as per the orthography profiles specification
            if line.lower().startswith("graphemes"):
                column_tokens = line.split("\t")

                # clean the header
                for column_token in column_tokens:
                    self.column_labels.append(column_token.lower().strip())
                continue

            # Unicode NFD the line
            line = unicodedata.normalize("NFD", line)

            # split the orthography profile into columns
            tokens = line.split("\t") 
            grapheme = tokens[0].strip()

            # check for duplicates in the orthography profile (fail if dups)
            if not grapheme in self.op_graphemes:
                self.op_graphemes[grapheme] = 1
            else:
                raise Exception("You have a duplicate in your orthography profile at: {0}".format(line_count))

            if len(tokens) == 1:
                continue

            # build transform table lookup
            for i in range(0, len(tokens)):
                token = tokens[i].strip()
                self.mappings[grapheme, self.column_labels[i].lower()] = token

            # move to rcParams
            # if rcParams['verbose']: print(rcParams['W_empty_cons'])
            # if rcParams['verbose']: print(rcParams['W_duplicates_in_op'])

        file.close()

        # print the trie structure if debug mode is on
        if self.debug:
            print("A graphical representation of your orthography profile in a trie ('*' denotes sentinels):\n")
            printTree(self.root, "")
            print()


    def _init_rules(self, file):
        """
        Process the orthography rules file.
        """
        try:
            codecs.open(file, 'r', 'utf-8')
        except IOError as e:
            print("\n[i] There is no file at the path you've specified.\n\n")

        rules_file = codecs.open(file, "r", 'utf-8')

        # compile the orthography rules
        for line in rules_file:
            line = line.strip()

            # skip any comments
            if line.startswith("#") or line == "":
                continue

            line = unicodedata.normalize("NFD", line)
            rule, replacement = line.split(",")
            rule = rule.strip() # just in case there's trailing whitespace
            replacement = replacement.strip() # because there's probably trailing whitespace!
            self.op_rules.append(re.compile(rule))
            self.op_replacements.append(replacement)
        rules_file.close()

        # check that num rules == num replacements; if not fail
        if len(self.op_rules) != len(self.op_replacements):
            raise ValueError("[i] Number of inputs does not match number of outputs in the rules file.")

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
        See: Unicode Standard Annex #29: UNICODE TEXT SEGMENTATION
        http://www.unicode.org/reports/tr29/

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
        """

        # init the regex Unicode grapheme cluster match
        grapheme_pattern = re.compile("\X", re.UNICODE)

        # add boundaries between words
        string = string.replace(" ", "#")

        # Unicode NDF the string
        string = unicodedata.normalize("NFD", string)

        result = ""
        graphemes = grapheme_pattern.findall(string)
        for grapheme in graphemes:
            result += grapheme+" "
        result = result.strip()
        return result


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
                if self.debug:
                    print("[i] The string '{0}' does not parse given the specified orthography profile.\n".format(word))
                    sys.stderr.write("The string '{0}' could not be tokenized".format(word))

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
        column = column.lower()

        if column == "graphemes":
            return self.graphemes(string)
        
        if not column in self.column_labels:
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

            # transform given the grapheme and column label
            result.append(self.mappings[token, column]) 

        return " ".join(result).strip()

    def transform_rules(self, string):
        """
        Convenience function that first tokenizes a string into orthographic profile-
        specified graphemes and then applies the orthography profile rules.
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

        # if no orthography profile rules file has been specified,
        # simply return the string
        if self.orthography_profile_rules == None:
            return string

        result = unicodedata.normalize("NFD", string)
        for i in range(0, len(self.op_rules)):
            match = self.op_rules[i].search(result)
            if not match == None:
                result = re.sub(self.op_rules[i], self.op_replacements[i], result)
                if self.debug:
                    print("input: ", string, "\t", "output: ", result)
                    print("pattern: ", self.op_rules[i].pattern, "\t", "replacement: ", self.op_replacements[i])
                    print()

        # this is incase someone introduces a non-NFD ordered sequence of characters
        # in the orthography profile
        result = unicodedata.normalize("NFD", result)
        return result


    def find_missing_characters(self, char_tokenized_string):
        """
        Given a string tokenized into characters, return a characters
        tokenized string where each character missing from the orthography 
        profile is replaced with a question mark <?>.
        """
        output = []
        chars = char_tokenized_string.split()
        for char in chars:
            if not char in self.op_graphemes:
                output.append("?")
            else:
                output.append(char)
        return " ".join(output)


    def combine_modifiers(self, string):
        """
        string == grapheme_clusters

        # Input should be tokenized on space...

        Given a string that is space-delimited on Unicode graphemes, group Unicode modifier letters with their preceeding base characters.

        Parameters
        ----------
        string : str
            A Unicode string to be parsed into Unicode cluster graphemes.

        .. todo:: check if we need to apply NDF after string is parsed

        """

        result = []
        graphemes = string[1:-1].split()
        temp = ""
        count = len(graphemes)

        for grapheme in reversed(graphemes):
            count -= 1
            if len(grapheme) == 1 and unicodedata.category(grapheme) == "Lm":
                temp = grapheme+temp

                # hack for the cases where a space modifier is the first character in the str
                if count == 0:
                    result[-1] = temp+result[-1]
                continue

            result.append(grapheme+temp)
            temp = ""
        return "# "+" ".join(result[::-1])+" #"

        
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

    def simple_ipa(self, string):
        """
        Simple IPA parsing
        """
        pass


# ---------- Tree node --------
class TreeNode(object):
    """
    Private class that creates the trie data structure from the orthography profile for parsing.
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

    file = codecs.open(file_name, "r", "utf-8")
    header = []

    for line in file:
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


