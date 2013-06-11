"""
Module to create sparse matrices and headers for Quantitative Language Comparison data.

"""

__author__ = "Steven Moran"
__date__ = "03-2013"

import sys
import collections
import gc

import numpy
from numpy import array
numpy.set_printoptions(threshold=numpy.nan) # set so everything will print

from scipy.sparse import csr_matrix, lil_matrix, coo_matrix, dia_matrix # do we need this to run this code?

import lingpy.sequence.ngram as ng # fix this ugly shit

from lingpy.basic.wordlist import Wordlist

class Matrix:
    """
    Module for creating sparse matrices from Lingpy wordlist data format.

    Parameters
    ----------

    input_data : str
        A string specifying the path and filename to the input data in QLC format.

    ngram_length : int (default = 2)
        The ngram window length; default is bigram.
        
        
    Notes
    -----

    The Matrix module 
    Takes as input a wordlist object that contains column data for language names, 
    concepts, and counterparts.
    """

    def __init__(self, filepath, ngram_length=2):
        """
        Matrix module to format wordlist data into various 2D matrices.
        """

        # TODO: add a debug parameter

        # Get a Wordlist object given the specified input file
        self.wl = Wordlist(filepath)
        self.ngram_length = ngram_length

        # check for language, concept, counterpart in Wordlist object; if missing data, fail
        self.wl_header = self.wl.header
        print(self.wl_header)

        # TODO: check for the items in the header
        """
        if all (k in self.wl_header for k in ("doculect", "concept", "orthoparse")):
            print("Matrix module input requires language ('doculect'), concept ('meaning'), and a qlc-format orthographic parse of the counterpart ('translation') and in wordlist object")
            sys.exit(1)

        sys.exit(1)
        """

        # if not, add one
        # print(self.wl.__getitem__(1))
        # print(self.wl[1,'orthoparse'])
        # a = lambda x:x.split('a')
        # wl.add_entries('juicyIPA','ipa',lambda x:x+x)
        
        # data structures to store various counts
        self._ngram_to_split_ngram = collections.defaultdict() # {"pb":"p_b"}

        # { word_id : { "#w_id":1, "wo_id":1, ...} } -- not ordered
        self._words_ngrams_counts = collections.defaultdict(lambda : collections.defaultdict(int)) 

        # { word_id : ["#w_id", "wo_id", ...] } -- ordered
        self._words_ngrams = collections.defaultdict(list)

        # data stored: {language: {counterpart: count} }
        self._languages_words_counts = collections.defaultdict(lambda : collections.defaultdict(int))

        # data stored: {concept: {counterpart: count} }
        self._concepts_words_counts = collections.defaultdict(lambda : collections.defaultdict(int))

        # data containers - using sets to discard duplicate ngrams
        non_unique_parsed_words = set()
        non_unique_ngrams = set()
        languages = set()
        concepts = set()
        unique_ngrams = set() 

        # loop over the wordlist data and parse into data strcutres
        for key in self.wl:
            language = self.wl[key, 'doculect']
            language = language.replace("_", "-") # fix
            concept = self.wl[key, 'concept']
            counterpart = self.wl[key,'orthoparse']
            # print(taxa, gloss, counterpart)

        # loop over the corpus reader data and parse into data structures

            
            """
            for language, concept, counterpart in language_concept_counterpart_iterator:
            # First do orthography parsing.
            if gram_type == "graphemes":
                parsed_counterpart_tuple = orthography_parser.parse_string_to_graphemes(counterpart) # graphemes
            elif gram_type == "phonemes":
                parsed_counterpart_tuple = orthography_parser.parse_string_to_ipa_phonemes(counterpart) # phonemes
            else:
                sys.exit('\ninvalid gram type: specify "phonemes" or "graphemes"\n')
                
            # TODO: move this to orthography parser
            # If string is unparsable, write to file.
            if parsed_counterpart_tuple[0] == False:
                invalid_parse_string = qlc.ngram.formatted_string_from_ngrams(parsed_counterpart_tuple[1])
                unparsables.write(language+"\t"+concept+"\t"+counterpart+"\t"+invalid_parse_string+"\n")
                continue
                """

            # parsed_counterpart = parsed_counterpart_tuple[1]

            counterpart = "# "+counterpart+" #"
            parsed_counterpart = tuple(counterpart.split())

            # Get ngrams as a tuple of tuples.
            # ngram_tuples = qlc.ngram.ngrams_from_graphemes(parsed_counterpart, ngram_length)
            ngram_tuples = ng.ngrams_from_graphemes(parsed_counterpart, ngram_length)

            # Format that tuple of tuples into a space-delimed string.
            ngrams_string = ng.formatted_string_from_ngrams(ngram_tuples)
            # print(ngrams_string)

            # Format tuple into unigrams split on "_" into a space-delimited string.
            split_ngrams_string = ng.split_formatted_string_from_ngrams(ngram_tuples)
            # print(split_ngrams_string)

            # check to make sure ngrams string ("#a ab b#") 
            # and split ngrams string ("#_a a_b b_#") are the same
            ngrams_string_list = ngrams_string.split()
            split_ngrams_string_list = split_ngrams_string.split()
            if len(ngrams_string_list) != len(split_ngrams_string_list):
                print("ngrams string and split ngrams sting do not match")
                sys.exit(1)

            # store key value pairs for ngram and split ngram; if unigram store the same
            for i in range(0, len(ngrams_string_list)):
                if self.ngram_length > 1:
                    self._ngram_to_split_ngram[ngrams_string_list[i]] = split_ngrams_string_list[i]
                else:
                    self._ngram_to_split_ngram[ngrams_string_list[i]] = ngrams_string_list[i]

            # Get the parsed version of counterparts.
            parsed_word = ng.formatted_string_from_ngrams(parsed_counterpart)
            # print("og: ", parsed_word)
            parsed_word = parsed_word.replace(" ", "")
            parsed_word = parsed_word.lstrip("#")
            parsed_word = parsed_word.rstrip("#")
            parsed_word = parsed_word.replace("#", " ")
            # print("pg: ", parsed_word)
            # print()

            # flipped
            # parsed_word_id = parsed_word+"_"+language
            parsed_word_id = language+"_"+parsed_word

            # if parsed_word not in dict:

            if not parsed_word_id in self._words_ngrams_counts:
                for ngram in ngrams_string.split():
                    # flipped
                    # non_unique_ngram = language+"_"+ngram
                    non_unique_ngram = language+"_"+ngram
                    non_unique_ngrams.add(non_unique_ngram)
                    self._words_ngrams_counts[parsed_word_id][non_unique_ngram] += 1
                    self._words_ngrams[parsed_word_id].append(non_unique_ngram)            


            # update data structures
            # self._languages_words_counts[language][parsed_word+"_"+language] += 1
            # self._concepts_words_counts[concept][parsed_word+"_"+language] += 1

            # flipped
            self._languages_words_counts[language][language+"_"+parsed_word] += 1
            self._concepts_words_counts[concept][language+"_"+parsed_word] += 1


            # add to header lists
            languages.add(language) # Append languages to unique set of langauge.
            concepts.add(concept)   # Append concepts to unique set of concepts.
            unique_ngrams.update(set(ngram_tuples)) # Append all the elements of ngram_tuples to unique_ngrams.


            # add to non-unique header lists
            # non_unique_parsed_words.add(parsed_word+"_"+language)
            # flipped
            non_unique_parsed_words.add(language+"_"+parsed_word)


        # listfy to sort
        self.languages = list(languages)
        self.languages.sort()

        self.concepts = list(concepts)
        self.concepts.sort()

        self.non_unique_parsed_words = list(non_unique_parsed_words)
        self.non_unique_parsed_words.sort()
        
        self.non_unique_ngrams = list(non_unique_ngrams)
        self.non_unique_ngrams.sort()

        self.unique_ngrams = list(unique_ngrams)
        self.unique_ngrams.sort()

    def get_wg_matrix(self):
        """
        Function returns a sparse matrix of language-specific words (rows) by language-specific graphemes (cols).
        The cells in the matrix contain the ngram count in the word.
        """
        wg = lil_matrix( (len(self.non_unique_parsed_words),len(self.non_unique_ngrams)), dtype=int)

        for i in range(0, len(self.non_unique_parsed_words)):
            # word_language_tokens = self.non_unique_parsed_words[i].partition("_")
            # word_language = word_language_tokens[2]
            # flipped
            language_word_tokens = self.non_unique_parsed_words[i].partition("_")
            language_word = language_word_tokens[2]

            for j in range(0, len(self.non_unique_ngrams)):
                # ngram_language_tokens = self.non_unique_ngrams[j].partition("_")
                # ngram_language = ngram_language_tokens[2]
                # flipped
                language_ngram_tokens = self.non_unique_ngrams[j].partition("_")
                language_ngram = language_ngram_tokens[2]

                # if the language IDs don't match, continue...
                # if word_language != ngram_language:
                #    continue
                # flipped
                if language_word != language_ngram:
                    continue

                # check to see if the word has the ngram
                if self._words_ngrams_counts[self.non_unique_parsed_words[i]][self.non_unique_ngrams[j]]:
                    wg[i,j] = self._words_ngrams_counts[self.non_unique_parsed_words[i]][self.non_unique_ngrams[j]]

                print("processing wg: "+str(i)+" "+str(j)+" of: "+str(len(self.non_unique_parsed_words))+" "+str(len(self.non_unique_ngrams)))
            # gc.collect()
        return wg

    def get_wg_coo_matrix(self):
        """
        Function returns a sparse matrix of language-specific words (rows) by language-specific graphemes (cols).
        The cells in the matrix contain the ngram count in the word.
        """
        rows = []
        cols = []
        data = []

        for i in range(0, len(self.non_unique_parsed_words)):
            rows.append(i)

            # word_language_tokens = self.non_unique_parsed_words[i].partition("_")
            # word_language = word_language_tokens[2]
            # flipped
            language_word_tokens = self.non_unique_parsed_words[i].partition("_")
            language_word = language_word_tokens[2]

            for j in range(0, len(self.non_unique_ngrams)):
                cols.append(j)
                # ngram_language_tokens = self.non_unique_ngrams[j].partition("_")
                # ngram_language = ngram_language_tokens[2]
                # flipped
                language_ngram_tokens = self.non_unique_ngrams[j].partition("_")
                language_ngram = language_ngram_tokens[2]

                # if the language IDs don't match, continue...
                # if word_language != ngram_language:
                #    continue
                # flipped
                if language_word != language_ngram:
                    continue



                # check to see if the word has the ngram
                if self._words_ngrams_counts[self.non_unique_parsed_words[i]][self.non_unique_ngrams[j]]:
                    data.append(self._words_ngrams_counts[self.non_unique_parsed_words[i]][self.non_unique_ngrams[j]])
                else:
                    data.append(0)
        # gc.collect()

        # wg = coo_matrix( (len(self.non_unique_parsed_words),len(self.non_unique_ngrams)), dtype=int)
        n_rows = array(rows)
        n_cols = array(cols)
        n_data = array(data)
        print(len(rows), len(cols), len(data))

        wg = coo_matrix((n_data, (n_rows,n_cols)), shape=( (len(self.non_unique_parsed_words),len(self.non_unique_ngrams) )))
        return wg

    # words/counterparts (rows) x languages (cols) x index (= if counterpart appears in that language)
    """
    def get_wl_matrix_by_language(self):
        # import scipy.linalg
        # a1 = np.array([[1,1,1],[1,1,1],[1,1,1]]
        # ...
        # scipy.linalg.block_diag(a1, a2, a3)

        #         ngram1_1 ngram2_1 ngram3_2 ngram4_2
        # word1_1
        # word2_1
        # word3_2
        # word4_2

        diagonal_matrices = []
        # block_diag((A, B, C)).todense()
        for language in self.languages:
            wl = lil_matrix(  (len(m._languages_words_counts[language]),1), dtype=int )

            for i in range(0, len(self.non_unique_parsed_words)):
                # word_language_tokens = self.non_unique_parsed_words[i].partition("_")
                # word_language = word_language_tokens[2]
                # flipped
                language_word_tokens = self.non_unique_parsed_words[i].partition("_")
                language_word = language_word_tokens[0]

                # print(language_word, language)
                for j in range(0, len(self.languages)):
                    # if the language IDs don't match, continue
                    if language_word != self.languages[j]:
                        continue

                    if self._languages_words_counts[self.languages[j]][self.non_unique_parsed_words[i]]:
                        wl[i,j] = 1
                gc.collect()
        return wl
        """

    # words/counterparts (rows) x languages (cols) x index (= if counterpart appears in that language)
    def get_wl_matrix(self):
        """
        Function returns a sparse matrix of words (rows) by languages (cols).
        The cells contain "1" if the word appears in the language.
        """
        wl = lil_matrix( (len(self.non_unique_parsed_words),len(self.languages)), dtype=int )

        for i in range(0, len(self.non_unique_parsed_words)):
            # Yorno-so_word
            # flipping
            # word_language_tokens = self.non_unique_parsed_words[i].partition("_")
            # word_language = word_language_tokens[2]

            language_word_tokens = self.non_unique_parsed_words[i].partition("_")
            language_word = language_word_tokens[2]

            language, separator, word = self.non_unique_parsed_words[i].partition("_")

            for j in range(0, len(self.languages)):
                # if the language IDs don't match, continue
                # flipping
                # if word_language != self.languages[j]:

                if language != self.languages[j]:
                    continue
                # print("w", language, "l", self.languages[j])

                if self._languages_words_counts[self.languages[j]][self.non_unique_parsed_words[i]]:
                    wl[i,j] = 1
            #gc.collect()
                print("processing wl: "+str(i)+" "+str(j)+" of: "+str(len(self.non_unique_parsed_words))+" "+str(len(self.languages)))

        return wl


    def get_wm_matrix(self):
        """
        Function returns a sparse matrix of words (rows) by meanings (aka concepts; cols).
        The cells contain "1" if the word appears with that meaning.
        
        """
        wm = lil_matrix( (len(self.non_unique_parsed_words),len(self.concepts)), dtype=int )

        for i in range(0, len(self.non_unique_parsed_words)):
            for j in range(0, len(self.concepts)):
                #if language_word != language_ngram:
                #    continue
                if self._concepts_words_counts[self.concepts[j]][self.non_unique_parsed_words[i]]:

                    # print(self._concepts_words_counts[self.concepts[j]][self.non_unique_parsed_words[i]])
                    # print(self.non_unique_parsed_words[i]+"\t"+self.concepts[j])

                    wm[i,j] = 1
            print("processing wm: "+str(i)+" "+str(j)+" of: "+str(len(self.non_unique_parsed_words))+" "+str(len(self.concepts)))
            # gc.collect()
        return wm

    def get_gp_matrix(self):
        gp = lil_matrix( (len(self.non_unique_ngrams),len(self.unique_ngrams)), dtype=int )
        for i in range(0, len(self.non_unique_ngrams)):
            for j in range(0, len(self.unique_ngrams)):
                # flipping
                # grapheme, separator, language_id = self.non_unique_ngrams[i].partition("_")
                language_id, separator, grapheme = self.non_unique_ngrams[i].partition("_")

                # get the phoneme-ngram; recompose from tuples
                phoneme = "" 
                for ngram in self.unique_ngrams[j]:
                    phoneme += ngram

                if grapheme == phoneme:
                    gp[i,j] = 1

            print("processing gp: "+str(i)+" "+str(j)+" of: "+str(len(self.non_unique_ngrams))+" "+str(len(self.unique_ngrams)))
            #gc.collect()
        return gp

    def write_header(self, list, source, ext):
        file = open(source+"/"+source+ext, "w")
        count = 0
        for item in list:
            count += 1
            file.write(str(count)+"\t"+item+"\n")
        file.close()

    def get_words_header(self):
        return self.get_header(self.non_unique_parsed_words)

    def get_meanings_header(self):
        return self.get_header(self.concepts)

    def get_ngrams_header(self):
        return self.get_header(self.non_unique_ngrams)

    def get_header(self, list):
        count = 0
        header = []
        for item in list:
            count += 1
            header.append(str(count)+"\t"+item)
        return header

    def get_split_ngrams_header(self, list):
        header = []
        for item in list:
            count, ngram_id = item.split("\t")
            # ngram, separator, id = ngram_id.partition("_")
            # flipped
            id, separator, ngram = ngram_id.partition("_")

            separated_ngram = self._ngram_to_split_ngram[ngram]
            # error check
            if ngram != separated_ngram.replace("_", ""):
                print("your grams aren't matching")
                sys.exit(1)
            header.append(count+"\t"+separated_ngram+"_"+id+"\t"+ngram_id)
        return header


    # makes phoneme_header ??
    def get_ngram_header(self):
        """
        Function to return a list of phonemes. 
        Count \t phoneme
        """
        count = 0
        header = []
        for i in range(0, len(self.unique_ngrams)):
            count += 1
            composed_ngram = ""
            for ngram in self.unique_ngrams[i]:
                composed_ngram += ngram
            header.append(str(count)+"\t"+composed_ngram)
        return header


    def get_words_ngrams_strings(self):
        return self._get_words_ngrams(False)

    def get_words_ngrams_indices(self):
        return self._get_words_ngrams(True)

    def _get_words_ngrams(self, index):
        header = []
        for word in self.non_unique_parsed_words:
            if not word in self._words_ngrams:
                print("warning: non_unique_parsed words does not match _words.ngrams")
                sys.exit(1)
            result = word

            for gram in self._words_ngrams[word]:
                # creates word-ngram indices: anene_7522 1 126 468 230 468 208
                if index:
                    gram = str(self.non_unique_ngrams.index(gram)+1) # add 1 because Python indexes from 0
                    result += "\t"+gram
                    continue

                # if not unigrams, delimit the graphemes on "_"
                if self.ngram_length > 1:
                    # creates word-ngram strings: anene_7522 #_a_7522 a_n_7522 n_e_7522 ... 
                    tokens = gram.partition("_")
                    # flipping
                    # unigram_item = self._ngram_to_split_ngram[tokens[0]]
                    unigram_item = self._ngram_to_split_ngram[tokens[2]]
                    # checks that the delimited ngrams are valid
                    # flipping
                    # if tokens[0] != unigram_item.replace("_", ""):
                    if tokens[2] != unigram_item.replace("_", ""):
                        print("your grams aren't matching")
                        sys.exit(1)
                    result += "\t"+unigram_item+tokens[1]+tokens[2]
                else:
                    result += "\t"+gram

            header.append(result)
        return header

if __name__=="__main__":
    m = Matrix("../../scripts/data/huber1992.csv")
