"""
This module provides graphemic and orthographic parsing with orthography profiles into QLC format.
"""

__author__ = "Steven Moran"
__date__ = "2010-12-01"

import sys
import unicodedata
import regex
import os

class DuplicateExceptation(Exception): pass

class GraphemeParser(object):
    """
    Class for Unicode graphemic parsing of Unicode strings. 

    Parameters
    ----------

    Notes
    -----
    This class handles (just) Unicode grapheme parsing, i.e. it uses the Unicode
    standard grapheme parsing rules, as implemented in the Python regex package 
    by Matthew Barnett, to do basic parsing with the "\X" grapheme regular 
    expression match. This grapheme match combines one or more Combining Diacritical 
    Marks to their base character. These are called "Grapheme clusters" in 
    Unicode parlance.

    This class is meant to do basic rudimentary parsing for things like getting 
    an additional unigram model (segments and their counts) in an input data source.

    For more elaborate orthographic parsing, see the OrthographyParser and the 
    OrthographyRulesParser.

    An additional method (in its infancy) called combine_modifiers handles the 
    case where there are Unicode Spacing Modifier Letters, which are not explicitly 
    combined to their base character in the Unicode Standard. These graphemes 
    are called "Tailored grapheme clusters" in Unicode. For more information 
    see the Unicode Standard Annex #29: Unicode Text Segmentation:

    http://www.unicode.org/reports/tr29/
    """

    def __init__(self):
        """
        Create the Unicode grapheme pattern match object.
        """
        self.grapheme_pattern = regex.compile("\X", regex.UNICODE)

    def combine_modifiers(self, string):
        """
        Given a string that is space-delimited on Unicode graphemes, group Unicode modifier letters with their preceeding base characters.

        Parameters
        ----------

        string : string
            A Unicode string to be parsed into graphemes.

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

    def parse_graphemes(self, string):
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
        """

        string = string.replace(" ", "#") # add boundaries between words
        string = unicodedata.normalize("NFD", string)
        result = "#"
        graphemes = self.grapheme_pattern.findall(string)
        for grapheme in graphemes:
            result += " "+grapheme
        result += " #"
        return (result)

    def parse_characters(self, string):
        """
        Given a string as input, return a space-delimited string of Unicode characters.

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
        result = "#"
        for character in string:
            result += " "+character
        result += " #"
        return (result)
        
    def parse_string_to_graphemes_string(self, string):
        """
        Deprecated function to parse graphemes and return a tuple of (T/F success, parse).

        .. todo:: Check code base for calls to this function.
        """
        string = string.replace(" ", "#") # add boundaries between words
        string = unicodedata.normalize("NFD", string)
        result = "#"
        graphemes = self.grapheme_pattern.findall(string)
        for grapheme in graphemes:
            result += " "+grapheme
        result += " #"
        # sys.stderr.write(result+"\n")
        return (True, result)

    def parse_string_to_graphemes(self, string):
        """
        Deprecated function to parse graphemes and return a tuple of (T/F success, tupled parse).

        .. todo:: Check code base for calls to this function.
        """
        (success, graphemes) = self.parse_string_to_graphemes_string(string)
        return (success, tuple(graphemes.split(" ")))

class OrthographyRulesParser(object):
    """
    Class for orthography rules parsing of Unicode strings.

    Parameters
    ----------
    orthography_profile_rules : file
        An orthography profile rules file.


    Notes
    -----
    Orthographic rules are current specified in a separate file from the 
    orthography profile (that specifics the document specific graphemes, 
    and possibly their IPA counterparts) and the orthography rules should 
    be applied to the output of an OrthographyParser.

    In an orthography rules file, rules are given in order in regular 
    expressions, e.g. this rule replaces a vowel followed by an <n> 
    followed by <space> followed by a second vowel with first vowel 
    <space> <n> <space> second vowel, e.g.:

    ([a|á|e|é|i|í|o|ó|u|ú])(n)(\s)([a|á|e|é|i|í|o|ó|u|ú]), \1 \2 \4

    .. todo:: integrate the rules file and the orthography profile into one.
    """

    def __init__(self, orthography_profile_rules):
        try:
            open(orthography_profile_rules)
        except IOError as e:
            print("\nWARNING: There is no file at the path you've specified.\n\n")

        self.rules = []
        self.replacements = []

        rules_file = open(orthography_profile_rules, "r", encoding="utf-8")

        # loop through the orthography fules and compile them
        for line in rules_file:
            line = line.strip()

            # skip any comments
            if line.startswith("#") or line == "":
                continue
            line = unicodedata.normalize("NFD", line)
            rule, replacement = line.split(",")
            rule = rule.strip() # just in case there's trailing whitespace
            replacement = replacement.strip() # because there's probably trailing whitespace!
            self.rules.append(regex.compile(rule))
            self.replacements.append(replacement)
        rules_file.close()

        # check that num rules == num replacements
        if len(self.rules) != len(self.replacements):
            print("there is a problem with your orthographic rules file: number of inputs does not match number of outputs")
            sys.exit(1)

    def parse_string(self, string):
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
        result = string
        for i in range(0, len(self.rules)):
            match = self.rules[i].search(result)
            if not match == None:
                result = regex.sub(self.rules[i], self.replacements[i], result)
        # this is incase someone introduces a non-NFD ordered sequence of characters
        # in the orthography profile
        result = unicodedata.normalize("NFD", result)
        return result

class OrthographyParser(object):
    """
    Class for orthographic parsing using orthography profiles as designed for the QLC project.

    Parameters
    ----------
    orthography_profile : file
        A document source-specific orthography profile.


    Notes
    -----
    The OrthographyParser reads in an orthography profile and calls a helper 
    class to build a trie data structure, which stores the possible Unicode 
    character combinations that are specified in the orthography profile 
    and appear in the data source.

    For example, an orthography profile might specify that in source X 
    <uu> is a single grapheme (Unicode parlance: tailored grapheme) and 
    thus should be chunked as so. Thus given an orthography profile and 
    some data to parse, the process would look like this:

    input string example: uubo
    output string example: # uu b o #

    where the output is given in QLC string format.

    Additionally, if a second column in an orthography profile is specified 
    (the first lists the graphemes in a given source), this class assumes 
    that that column is the IPA translation of the graphemes. A dictionary 
    is created that keeps a mapping between source-specific graphemes and 
    their IPA counterparts.

    Deprecated methods in this class return a tuple of (True or False, parsed-string).
    The first element in the tuple relays whether the string parsed sucessfully. The 
    second element returns the parsed string.

    .. todo:: remove any code that uses the deprecated methods
    .. todo:: update the processes so we don't just assume col 2 is IPA
    """

    def __init__(self, orthography_profile):
        try:
            open(orthography_profile)
        except IOError as e:
            print("\nWARNING: There is no file at the path you've specified.\n\n")

        # path and filename of orthography profile
        self.orthography_profile = orthography_profile

        # read in orthography profile and create a tree structure
        self.root = createTree(orthography_profile)

        # lookup table
        self.grapheme_to_phoneme = {}

        # create look up table of grapheme to IPA from orthography profile
        file = open(orthography_profile, "r", encoding="utf-8")

        # an orthography profile may have one or more comma delimited columns
        # the first column is graphemes, the second IPA, etc.
        self.multiple_columns = False

        line_count = 0
        for line in file:
            line_count += 1
            line = line.strip()

            # skip any comments
            if line.startswith("#") or line == "":
                continue

            line = unicodedata.normalize("NFD", line)

            tokens = line.split(",") # split the orthography profile into columns

            if len(tokens) > 1:
                self.multiple_columns = True

            grapheme = tokens[0].strip()
            phoneme = tokens[1].strip()

            # TODO: assuming that we have two filled columns; we should throw an error otherwise
            if not grapheme in self.grapheme_to_phoneme:
                self.grapheme_to_phoneme[grapheme] = phoneme
            else:
                raise DuplicateException("You have a duplicate in your orthography profile at: {0}".format(line_count))
        file.close()

        # uncomment this line if you want to see the orthography profile tree structure
        # printTree(self.root, "")

    def exists_multiple_columns(self):
        """
        Returns boolean of whether multiple columns exist in the orthography profile, e.g. graphemes and IPA and x, etc.
        """
        return self.multiple_columns

    def parse_string_to_graphemes_string_DEPRECATED(self, string):
        """
        Deprecated function to parse str into tuples (success, parsed str).
        """
        string = string.replace(" ", "#") # add boundaries between words
        string = unicodedata.normalize("NFD", string)
        result = ""
        result += printMultigraphs(self.root, string, result+"# ")
        return (True, result)

    def parse_graphemes(self, string):
        """
        Parses orthograhy profile specified graphemes given a string.

        Parameters
        ----------
        string : str
            The str to be parsed and formatted.

        Returns
        -------
        result : str
            The result of the parsed and QLC formatted str.

        Notes
        -----
        An example: 
        
        input: dog shit
        output # d o g # sh i t #

        .. todo:: updated for LINGPY: deal with error handling in THIS class
        .. todo:: write unparsables to disk!
        """

        success = True
        parses = []
        string = unicodedata.normalize("NFD", string)
        for word in string.split():
            parse = getParse(self.root, word)
            if len(parse) == 0:
                print("[i] The string <"+string+"> does not parse given the specified orthography profile.\n")
                parse = " <no-valid-parse> "
            parses.append(parse)
        # Use "#" as a word boundary token (a special 'grapheme').
        result = "".join(parses).replace("##", "#")  # Ugly. Oh well.
        return result


    # updated for LINGPY - deal with error handling in THIS class
    # TODO -- write unparsables to disk!
    def graphemes_to_ipa(self, string):
        """
        Returns the parsed and formated string given the orthography profile.

        Parameters
        ----------
        string : str
            The str to be parsed and formatted.

        Returns
        -------
        ipa : str
            The str parsed and formatted and flipped into IPA.

        Notes
        -----
        Graphemes encoded in the orthography profile and the IPA row are used 
        as a global scope lookup using a dict.
        """
        graphemes = self.parse_graphemes(string)
        ipa = graphemes

        # flip the graphemes into phonemes
        # TODO: probably don't need a loop for *every string* -- refactor
        # ASSUMES that every grapheme has a corresponding ipa character

        for k, v in self.grapheme_to_phoneme.items():
            if k == "#" or k == " ":
                continue
            ipa = ipa.replace(k, v)

        return ipa

    def parse_string_to_graphemes_string(self, string):
        """
        Deprecated methods that returns parsed str in a tuple.

        .. todo:: check for code calling this function then remove
        """
        success = True
        parses = []
        string = unicodedata.normalize("NFD", string)
        for word in string.split():
            # print("word: "+"\t"+word)
            parse = getParse(self.root, word)
            if len(parse) == 0:
                success = False
                # parse = "# <no valid parse> #"
                # parse = word
                parse = " <no-valid-parse> "

            parses.append(parse)

        # Use "#" as a word boundary token (a special 'grapheme').
        result = "".join(parses).replace("##", "#")  # Ugly. Oh well.
        return (success, result)

    def parse_string_to_graphemes(self, string):
        """
        Deprecated function that parses string and returns tuple of graphemes.

        .. todo:: check code calling this function and remove
        """
        (success, graphemes) = self.parse_string_to_graphemes_string(string)

        return (success, tuple(graphemes.split(" ")))

    def parse_string_to_ipa_phonemes(self, string):
        """
        Deprecated function to parse string and returns tuple of success and phonemes.

        .. todo:: check for code calling this funtion and update and remove
        """
        (success, graphemes) = self.parse_string_to_graphemes_string(string)
        if not success:
            return (False, graphemes)

        # flip the graphemes into phonemes
        # this is so ghetto and fragile -- depends on the precise encoding of the orthography profile
        
        graphemes = graphemes.split(" ")
        ipa = []
        for i in range (0, len(graphemes)):
            if graphemes[i] == "#":
                ipa.append(graphemes[i])
                continue
            grapheme = self.grapheme_to_phoneme[graphemes[i]]
            if grapheme != "" and grapheme != " ":
                ipa.append(grapheme)

        return (success, tuple(ipa))

    def parse_string_to_ipa_string(self, string):
        """
        Deprecated function to parse str into tuple of success and phonemes.

        .. todo:: check for caller code and update
        """
        (success, graphemes) = self.parse_string_to_graphemes_string(string)
        if not success:
            return (False, graphemes)
        ipa = graphemes

        # flip the graphemes into phonemes
        # TODO: probably don't need a loop for *every string* -- refactor
        for k, v in self.grapheme_to_phoneme.items():
            ipa = ipa.replace(k, v)

        return (True, ipa)

    def parse_formatted_string_to_ipa_string(self, string):
        """
        Deprecated function to parse formatted string into graphemes.

        .. todo:: check for caller code and update
        """
        string = string.strip()
        # flip the graphemes into phonemes
        # TODO: probably don't need a loop for *every string* -- refactor
        # TODO: this method assumes exceptions have already been caught by
        #  orthographic and orthographic rule parsings
        result = ""
        for char in string.split():
            if char == "#":
                result += " "+char
                continue
            if not char in self.grapheme_to_phoneme:
                print("you are missing the grapheme in your orthography profile!")
                sys.exit(1)
            if not self.grapheme_to_phoneme.__contains__(char):
                print("you are missing the grapheme in your orthography profile!!!")
                sys.exit(1)
            result += " "+self.grapheme_to_phoneme[char]
        return result
    

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

    #file = codecs.open(file_name, "r", "utf-8")
    file = open(file_name, "r", encoding="utf-8")

    for line in file:
        #print(line.encode("utf-8"))
        line = line.strip()

        # skip any comments
        if line.startswith("#") or line == "":
            continue

        line = unicodedata.normalize("NFD", line)
        tokens = line.split(",") # split the orthography profile into columns
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

# ---------- Main ------
if __name__=="__main__":
    o = OrthographyParser("../data/orthography_profiles/thiesen1998.txt")
    g = GraphemeParser()
    test_words = ["aa", "aabuu", "uuabaa auubaa"]
    print()
    for word in test_words:
        print("original word:", word)
        print("parse_string_to_graphemes_string:", o.parse_string_to_graphemes_string(word))
        print("parse_string_to_ipa_string:", o.parse_string_to_ipa_string(word))
        print("parse_string_to_graphemes:", o.parse_string_to_graphemes(word))
        # print("parse_graphemes:", g.parse_graphemes(word))
        print()

