# -*- coding: utf-8 -*-
"""
Module to create sparse matrices and headers for Quantitative Language Comparison data.

"""
import sys
import collections
import gc
import qlc.ngram
from qlc.corpusreader import CorpusReaderWordlist
import numpy
from numpy import array
from scipy.sparse import csr_matrix, lil_matrix, coo_matrix, dia_matrix # do we need this to run this code?
numpy.set_printoptions(threshold=numpy.nan) # set so everything will print

class Matrix:
    """
    Basic class for creating matrices for QLC data.
    """

    def __init__(self, language_concept_counterpart_iterator, orthography_parser, gram_type, ngram_length):

        self.ngram_length = ngram_length
        
        # write to disk the forms that can't be parsed
        unparsables = open("output/unparsables.txt", "w")

        # data structures
        self._ngram_to_split_ngram = collections.defaultdict() # {"pb":"p_b"}
        
        # { word_id : { "#w_id":1, "wo_id":1, ...} } -- not ordered
        self._words_ngrams_counts = collections.defaultdict(lambda : collections.defaultdict(int)) 

        # { word_id : ["#w_id", "wo_id", ...] } -- ordered
        self._words_ngrams = collections.defaultdict(list)

        # data stored: {language: {counterpart: count} }
        self._languages_words_counts = collections.defaultdict(lambda : collections.defaultdict(int))

        # data stored: {concept: {counterpart: count} }
        self._concepts_words_counts = collections.defaultdict(lambda : collections.defaultdict(int))

        # data stored: {language: {concept:1} }
        self._language_concepts = collections.defaultdict(lambda : collections.defaultdict(int))



        # data containers
        non_unique_parsed_words = set()
        non_unique_ngrams = set()
        languages = set()
        concepts = set()
        unique_ngrams = set() # using a Python set discards duplicate ngrams

        # loop over the corpus reader data and parse into data structures
        for language, concept, counterpart in language_concept_counterpart_iterator:
            # First do orthography parsing.
            if gram_type == "graphemes":
                parsed_counterpart_tuple = orthography_parser.parse_string_to_graphemes(counterpart) # graphemes
            elif gram_type == "phonemes":
                parsed_counterpart_tuple = orthography_parser.parse_string_to_ipa_phonemes(counterpart) # phonemes
            elif gram_type == "already_parsed":
                grams_list = counterpart.split()
                grams_tuple = tuple(grams_list)
                parsed_counterpart_tuple = (True, grams_tuple)
                
            else:
                sys.exit('\ninvalid gram type: specify "phonemes" or "graphemes"\n')
                
            # TODO: move this to orthography parser
            # If string is unparsable, write to file.
            if parsed_counterpart_tuple[0] == False:
                invalid_parse_string = qlc.ngram.formatted_string_from_ngrams(parsed_counterpart_tuple[1])
                unparsables.write(language+"\t"+concept+"\t"+counterpart+"\t"+invalid_parse_string+"\n")
                continue
            parsed_counterpart = parsed_counterpart_tuple[1]

            # Get ngrams as a tuple of tuples.
            ngram_tuples = qlc.ngram.ngrams_from_graphemes(parsed_counterpart, ngram_length)

            # Format that tuple of tuples into a space-delimed string.
            ngrams_string = qlc.ngram.formatted_string_from_ngrams(ngram_tuples)
            # print(ngrams_string)

            # Format tuple into unigrams split on "_" into a space-delimited string.
            split_ngrams_string = qlc.ngram.split_formatted_string_from_ngrams(ngram_tuples)
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
            parsed_word = qlc.ngram.formatted_string_from_ngrams(parsed_counterpart)
            # print("og: ", parsed_word)
            parsed_word = parsed_word.replace(" ", "")
            parsed_word = parsed_word.lstrip("#")
            parsed_word = parsed_word.rstrip("#")
            parsed_word = parsed_word.replace("#", " ")
            # print("pg: ", parsed_word)
            # print()

            parsed_word_id = parsed_word+"_"+language

            # if parsed_word not in dict:

            if not parsed_word_id in self._words_ngrams_counts:
                for ngram in ngrams_string.split():
                    non_unique_ngram = language+"_"+ngram
                    non_unique_ngrams.add(non_unique_ngram)
                    self._words_ngrams_counts[parsed_word_id][non_unique_ngram] += 1
                    self._words_ngrams[parsed_word_id].append(non_unique_ngram)            


            # update data structures
            self._languages_words_counts[language][parsed_word+"_"+language] += 1
            self._concepts_words_counts[concept][parsed_word+"_"+language] += 1
            self._language_concepts[language][concept] += 1

            # add to header lists
            languages.add(language) # Append languages to unique set of langauge.
            concepts.add(concept)   # Append concepts to unique set of concepts.
            unique_ngrams.update(set(ngram_tuples)) # Append all the elements of ngram_tuples to unique_ngrams.


            # add to non-unique header lists
            non_unique_parsed_words.add(parsed_word+"_"+language)


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
            word, separator, language = self.non_unique_parsed_words[i].partition("_")
            # word_language_tokens = self.non_unique_parsed_words[i].partition("_")
            # language
            word_language = language

            for j in range(0, len(self.non_unique_ngrams)):
                # ngram_language_tokens = self.non_unique_ngrams[j].partition("_")
                language, separator, ngram = self.non_unique_ngrams[j].partition("_")
                # print(ngram_language_tokens)
                # ngram_language = ngram_language_tokens[2]
                ngram_language = language

                # if the language IDs don't match, continue...
                if word_language != ngram_language:
                    # print(word_language)
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
            word_language_tokens = self.non_unique_parsed_words[i].partition("_")
            word_language = word_language_tokens[2]

            for j in range(0, len(self.non_unique_ngrams)):
                cols.append(j)
                ngram_language_tokens = self.non_unique_ngrams[j].partition("_")
                ngram_language = ngram_language_tokens[2]

                # if the language IDs don't match, continue...
                if word_language != ngram_language:
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
            word_language_tokens = self.non_unique_parsed_words[i].partition("_")
            word_language = word_language_tokens[2]
            word, separator, language = self.non_unique_parsed_words[i].partition("_")

            for j in range(0, len(self.languages)):
                # if the language IDs don't match, continue
                if word_language != self.languages[j]:
                    continue
                # here
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
            word, separator, language = self.non_unique_parsed_words[i].partition("_")
            for j in range(0, len(self.concepts)):
                if not self.concepts[j] in self._language_concepts[language]:
                    continue
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
                language_id, separator, grapheme = self.non_unique_ngrams[i].partition("_")
                # print(grapheme, separator, language_id)
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
                    language, separator, ngram = gram.partition("_")
                    unigram_item = self._ngram_to_split_ngram[ngram]

                    # checks that the delimited ngrams are valid
                    if ngram != unigram_item.replace("_", ""):
                        print("your grams aren't matching")
                        sys.exit(1)
                    result += "\t"+unigram_item+separator+language
                else:
                    result += "\t"+gram

            header.append(result)
        return header







if __name__=="__main__":
    import sys
    import time
    from qlc.corpusreader import CorpusReaderWordlist
    from qlc.orthography import OrthographyParser, GraphemeParser
    from scipy.io import mmread, mmwrite # write sparse matrices

    if len(sys.argv) != 3:
        print("call: python matrix.py bibtex_key output_directory\n")
        print("python matrix.py huber1992 huber1992-bigrams\n")

    bibtex_key = sys.argv[1] # bibtex key of source in the database, e.g. zgraggen1980
    source = sys.argv[2] # dictionary/wordlist source key
    output_dir = "output/"+source+"/"

    # get data from corpus reader
    cr = CorpusReaderWordlist("data/csv")          # real data
    # cr = CorpusReaderWordlist("data/testcorpus") # test data

    # initialize orthography parser for source
    o = OrthographyParser("data/orthography_profiles/"+source+".txt")
    # o = GraphemeParser() # or use the grapheme parser

    # create a generator of corpus reader data
    wordlist_iterator = ( (wordlistdata_id, concept, counterpart)
        for wordlistdata_id in cr.wordlistdata_ids_for_bibtex_key(bibtex_key)
        for concept, counterpart in cr.concepts_with_counterparts_for_wordlistdata_id(wordlistdata_id)
    )

    # write the data to disk -- note it exhausts the generator, so either the generator
    # must be "regenerated" or run the following lines without the rest of the code below
    # move this into a method in the class

    """
    file = open(output_dir+source+"_data.txt", "w")    
    file.write("# wordlistdata_id"+"\t"+"language bookname"+"\t"+"concept"+"\t"+"counterpart"+"\n")
    for wordlistdata_id, concept, counterpart in wordlist_iterator:
        result = wordlistdata_id+"\t"+cr.get_language_bookname_for_wordlistdata_id(wordlistdata_id)+"\t"+concept+"\t"+counterpart+"\n"
        file.write(result)
    file.close()
    """

    # initialize matrix class
    m = Matrix(wordlist_iterator, o, "phonemes", 2) # pass ortho parser and ngram length

    t0=time.time()

    # write the sparse matrices to disk...

    wg = m.get_wg_matrix()
    mmwrite(output_dir+source+"_WG.mtx", wg)

    wl = m.get_wl_matrix()
    mmwrite(output_dir+source+"_WL.mtx", wl)

    wm = m.get_wm_matrix()
    mmwrite(output_dir+source+"_WM.mtx", wm)

    gp = m.get_gp_matrix()    
    mmwrite(output_dir+source+"_GP.mtx", gp)


    # write header files...

    # 1. write header: wordlistid \t language name
    file = open(output_dir+source+"_wordlistids_lgnames_header.txt", "w")
    wordlist_ids = []
    for wordlistdata_id in cr.wordlistdata_ids_for_bibtex_key(bibtex_key):
        wordlist_ids.append(wordlistdata_id)
        # print(len(wordlist_ids))
        
    wordlist_ids.sort()
    for wordlist_id in wordlist_ids:
        file.write(wordlist_id+"\t"+cr.get_language_bookname_for_wordlistdata_id(wordlist_id)+"\n")
    file.close()


    # 2. write words header: count \t word_id 
    file = open(output_dir+source+"_words_header.txt", "w")    
    words_header = m.get_words_header()
    for item in words_header:
        file.write(item+"\n")
    file.close()
    
    # 3. write meanings header: count \t meaning (concept)
    file = open(output_dir+source+"_meanings_header.txt", "w")    
    meanings_header = m.get_meanings_header()
    for item in meanings_header:
        file.write(item+"\n")
    file.close()

    # 4. write ngrams header: count \t ngram
    file = open(output_dir+source+"_ngrams_header.txt", "w")    
    ngrams_header = m.get_ngrams_header()
    # additional call to add in "_" separated ngrams
    split_ngrams_header = m.get_split_ngrams_header(ngrams_header)
    for item in split_ngrams_header:
        file.write(item+"\n")
    file.close()

    # 5. write phonemes header: count \t phoneme
    file = open(output_dir+source+"_phonemes_header.txt", "w")        
    phoneme_header = m.get_ngram_header()
    for item in phoneme_header:
        file.write(item+"\n")
    file.close()

    # 6. write the word and ngrams/ngrams-indices
    file = open(output_dir+source+"_words_ngrams_strings.txt", "w")
    words_ngrams_strings_header = m.get_words_ngrams_strings()    
    for item in words_ngrams_strings_header:
        file.write(item+"\n")
    file.close()

    file = open(output_dir+source+"_words_ngrams_indices.txt", "w")
    words_ngrams_indices_header = m.get_words_ngrams_indices()
    for item in words_ngrams_indices_header:
        file.write(item+"\n")
    file.close()

    print("total time to fill the matrices:", time.time() -t0)
