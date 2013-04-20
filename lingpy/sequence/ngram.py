"""
This module provides some basic ngram functions, such as parsing QLC-string formatted strings into ngrams sequences and generating unigram models for intial orthography profiles.
"""

__author__ = "Steven Moran"
__date__ = "2011-01-01"

import collections
import numpy
import sys
import operator
import unicodedata

# don't use absolute names for lingpy-imports
from .orthography import OrthographyParser

def character_model(list):
    """
    Method counts characters and their frequences given a list of words.
    """
    segments_hash = collections.defaultdict(int)
    segment_count = 0

    for word in list:
        word = word.strip()

        for char in word:
            segment_count += 1
            segments_hash[segment] += 1            

    # sort the hash set by values (frequencies of unigrams)
    segments_sorted = sorted(segments_hash.items(), key=operator.itemgetter(1), reverse=True)
    
    print("Character"+"\t"+"Count"+"\t"+"Frequency") 
    for segment in segments_sorted:
        segment, count = segment[0], segment[1]
        frequency = segments_hash[segment]/segment_count
        print(segment+"\t"+str(count)+"\t"+str(frequency))

def unigram_model(list):
    """
    Takes as input a list of QLC tokenized words and outputs a unigram model.
    """
    segments_hash = collections.defaultdict(int)
    segment_count = 0

    for word in list:
        word = word.strip()
        segments = word.split()
        for segment in segments:
            segment_count += 1
            segments_hash[segment] += 1

    # sort the hash set by values (frequencies of unigrams)
    segments_sorted = sorted(segments_hash.items(), key=operator.itemgetter(1), reverse=True)
    
    print("Phone"+"\t"+"Count"+"\t"+"Frequency") # +"\t"+"plog")
    for segment in segments_sorted:
        segment, count = segment[0], segment[1]
        frequency = segments_hash[segment]/segment_count
        plog = "coming..."
        print(segment+"\t"+str(count)+"\t"+str(frequency)) # +"\t"+str(plog))

def unicode_categories(list):
    """
    Return the type of unicode categories and their characters, e.g.

    Ll (Letter, Lowercase): a, b, c, d, e...
    Lm (Letter, Modifier): ', :, ...
    """
    pass

def extended_unicode_model(list):
    """
    Takes as input a list of QLC-formatted words and outputs a unigram model.
    """
    segments_hash = collections.defaultdict(int)
    segment_count = 0

    for word in list:
        word = word.strip()
        segments = word.split()
        for segment in segments:
            segment_count += 1
            segments_hash[segment] += 1

    segments_sorted = sorted(segments_hash.items(), key=operator.itemgetter(1), reverse=True)
    
    # print("Phone"+"\t"+"Int"+"\t"+"Count"+"\t"+"Frequency") # +"\t"+"plog")
    print("Char"+"\t"+"int"+"\t"+"Unicode name"+"\t"+"category"+"\t"+"comb class"+"\t"+"decomposition"+"\t"+"count"+"\t"+"frequency")

    for segment in segments_sorted:
        segment, count = segment[0], segment[1]
        frequency = segments_hash[segment]/segment_count

        # decimal = unicodedata.decimal(segment)
        name = unicodedata.name(segment)
        category = unicodedata.category(segment)
        combining_class = unicodedata.combining(segment)
        decomposition = unicodedata.decomposition(segment)

        print(segment+"\t"+str(ord(segment))+"\t"+name+"\t"+category+"\t"+str(combining_class)+"\t"+decomposition+"\t"+str(count)+"\t"+str(frequency))


def ngrams_from_graphemes(graphemes, n=1):
    """
    Takes a tuple of (orthographically parsed) graphemes and returns a tuple of ngrams.

    Parameters
    ----------
    graphemes : tuple
        a tuple of grapheme strings from which the tuple of ngrams is extracted

    n : int (default = 1)
        the number of graphemes that have to be looked at for the ngram
        
    Return
    ------
    _ : tuple
        a tuple of tupled ngrams for the input string

    Notes
    -----
    Default is 1 ngram, i.e. methods will return a tuple of unigrams.

    Note: a word in a dataset should never have less than three elements, 
    e.g. if a word is composed of only one sound (represented by one grapheme, 
    then with the word boundaries ("#"), the length of the word should 
    be at least 3: "#u#". Therefore a tuple of length three will be returned: 
    ("#", "u", "#"). This script will fail if the input tuple is length less 
    than 3.

    Note also that one element tuples (as casted from a list) contain a comma, e.g.:
    >>> l = ["1"]
    >>> print tuple(l)
    ('1',)
    """

    # we don't except ngrams less than length 3 (1 sound + 2 word boundaries)
    if len(graphemes) < 3:
        print()
        print("WARNING: the qlc.ngram class says that you have a word (inclusive of word boundaries '#') that is less than length 3. That's bad!")
        print("word: "+graphemes+"\nExiting...\n")
        sys.exit(1)

    list_of_grams = []

    # if more than unigrams
    if n > 1:
        for i in range(0, len(graphemes)-n+1):
            # print("graphemes[i:i+n]", graphemes[i:i+n])
            list_of_grams.append(graphemes[i:i+n])    
            # print("list_of_grams:", list_of_grams)
            # print()
    else:
        for i in range(0, len(graphemes)):
            list_of_grams.append(tuple(graphemes[i]))
    return tuple(list_of_grams)

"""
    print(type(graphemes))
    list_of_grams = []
    for i in range(0, len(graphemes)-n+1):
        print(graphemes[i:i+n])
        list_of_grams.append(graphemes[i:i+n])
    return list_of_grams
    """

def formatted_string_from_ngrams(ngrams_tuple):
    """ Convert a tuple of tuples to a list of strings, and return joined string. """
    ngrams_list = []
    for ngram in ngrams_tuple:
        ngrams_list.append("".join(list(ngram)))
    return " ".join(ngrams_list)

def split_formatted_string_from_ngrams(ngrams_tuple):
    """ Convert a tuple of tuples to a list of strings, split on unigrams, and return joined string. """
    ngrams_list = []
    for ngram in ngrams_tuple:
        ngrams_list.append("_".join(list(ngram)))
    return " ".join(ngrams_list)

def words_ngrams_list_for_graphemes_list(graphemes_list, n=1):
    ngrams_list = []
    for graphemes in graphemes_list:
        ngrams = ngrams_from_graphemes(graphemes, n)
        ngrams_list.extend(ngrams)
    return ngrams_list

def words_ngrams_matrix_for_graphemes_list(graphemes_list, n=1):
    """
    Goes through the list of graphemes and gathers all the ngrams
    
    Parameters
    ----------
    grapheme_list: list of strings
        a list of strings from which the list of ngrams is extracted
    n: integer
        the number of graphemes that have to be looked at for the ngram
        default: n = 2 (bigram mode)
        
    Returns
    -------
    row_names: list
        list of strings from the input list (row names of the matrix)
    ngrams_list: list
        list of ngrams that occur in the matrix
    matrix: array
        a matrix of ngrams with the list of graphemes as its rows, the
        list of n-grams as its columns and the counts of each ngram for
        the respective word in the cells
    """
    ngrams_counts = list()
    ngrams_set = set()
    row_names = list()
    # go through the list of words
    for i, graphemes in enumerate(graphemes_list):

        ngrams_counts.append(collections.defaultdict(int))
        row_names.append(graphemes)
        ngrams = ngrams_from_graphemes(graphemes, n)

        # go through the list of ngrams and store them in a set
        for ngram in ngrams:
            ngrams_counts[i][ngram] += 1
            ngrams_set.add(ngram)

    ngrams_list = sorted(list(ngrams_set))
    # generate a matrix with zeros

    matrix = numpy.zeros( (len(row_names), len(ngrams_list)) )
    # matrix = numpy.zeros( ( len(row_names), len(ngrams_list) ) )

    
    # fill the matrix with the ngram counts
    for i in range(len(row_names)):
        for j, ngram in enumerate(ngrams_list):
            matrix[i][j] = ngrams_counts[i][ngram]

    return matrix
        
if __name__ == '__main__':
    # NgramTest().run()
    # tuple = ('#', 'h', 'a', 'd', 'ɯ', '#')
    graphemes = ('#', 'h', '#')
    # string = "#h#"
    print()
    print("ngrams_from_graphemes, tuple:", ngrams_from_graphemes(graphemes, 1))
    print()

#    from lingpy.sequence import tokenizer
#    t = tokenizer.tokenizer()
#    unigram_model(t)
