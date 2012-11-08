================================
Handling Multilingual Word Lists
================================

What is a Word List?
--------------------

Generally, a word list is a simple tabular data structure in which multiple
languages are structured in such a way that **words** are ordered in rows and
columns according to the **language** to which they belong and the **concept** they
denote. A simple word list could thus be displayed as a simple tab-delimited
text file::

    CONCEPT     GERMAN      ENGLISH     RUSSIAN     UKRAINIAN
    hand        Hand        hand        рука        рука
    leg         Bein        leg         нога        нога
    Woldemort   Waldemar    Woldemort   Владимир    Володимир
    Harry       Harald      Harry       Гарри       Гаррi

However, this format has a striking drawback, in so far as what we call a "word"
can have multiple manifestations in our data. Thus, the same word list could
look like this, if we preferred to have the words represented in phonetic
transcription::

    CONCEPT     GERMAN      ENGLISH     RUSSIAN     UKRAINIAN
    hand        hant        hænd        ruka        ruka
    leg         bain        lɛg         noga        noga
    Woldemort   valdəmar    wɔldəmɔrt   vladimir    volodimir
    Harry       haralt      hæri        gari        gari

And if we wanted to display only which of the words are cognate, we could
represent it in a numerical format, where all words sharing the same number are
cognate::

    CONCEPT     GERMAN      ENGLISH     RUSSIAN     UKRAINIAN
    hand        1           1           2           2
    leg         3           4           5           5  
    Woldemort   6           6           6           6
    Harry       7           7           7           7

When dealing with word lists in general, we thus need something more than just a
two-dimensional representation format. A solution is to use a simple csv-format
with a header which specifies not only the concept and the language, but also
all different possible **entry-types** a word can have::

    @author: Potter, Harry
    @date: 2012-11-07
    ID   CONCEPT     COUNTERPART   IPA         DOCULECT     COGID
    1    hand        Hand          hant        German       1
    2    hand        hand          hænd        English      1
    3    hand        рука          ruka        Russian      2
    4    hand        рука          ruka        Ukrainian    2
    5    leg         Bein          bain        German       3
    6    leg         leg           lɛg         English      4
    7    leg         нога          noga        Russian      5
    8    leg         нога          noga        Ukrainian    5
    9    Woldemort   Waldemar      valdemar    German       6
    10   Woldemort   Woldemort     wɔldemɔrt   English      6
    11   Woldemort   Владимир      vladimir    Russian      6
    12   Woldemort   Володимир     volodimir   Ukrainian    6
    9    Harry       Harald        haralt      German       7
    10   Harry       Harry         hæri        English      7
    11   Harry       Гарри         gari        Russian      7
    12   Harry       Гаррi         gari        Ukrainian    7

This format is, of course, much more redundant, than the word list format, but
it allows to display multiple entry-types for the counterparts of a given
concept in a given language. Moreover, this format is the basic of the
py:class::`~lingpy.basic.wordlist.WordList` class in LingPy, which makes it easy
to handle word lists with multiple entry-types of words.

The above-given csv-file `harry_potter.csv` is available in the test folder of LingPy.
In order to get it loaded, we first call load it into a dictionary::

    >>> d = _load_dict('harry_potter.csv')

Then we can use this dictionary to create our WordList instance::

    >>> wl = WordList(d)

If one wants to access only the IPA values in tabular format, all one has to do
is::

    >>> wl.ipa
    [['wɔldemɔrt', 'valdemar', 'vladimir', 'volodimir'],
     ['hæri', 'haralt', 'gari', 'gari'],
     ['lɛg', 'bain', 'noga', 'noga'],
     ['hænd', 'hant', 'ruka', 'ruka']]

The same for cognates::

    >>> wl.cognate
    [[6, 6, 6, 6], [7, 7, 7, 7], [4, 3, 5, 5], [1, 1, 2, 2]]

Or for the languages and the concepts in the dataset::

    >>> wl.language
    ['English', 'German', 'Russian', 'Ukrainian']
    >>> wl.concept
    ['Harry', 'Woldemort', 'hand', 'leg']
    


    


How are Word Lists defined?
---------------------------

In LingPy the WordList class handles wordlists. 

The wordlist.conf file
----------------------

The structure of word lists is defined by the configuration file `wordlists.rc`. This file is
automatically loaded when initializing a WordList instance:

    >>> wl = WordList(data)

It can, however, also be passed by the user:

    >>> wl = WordList(data,conf="path_to_file")

The file is a simple tab-delimited csv-file and has the following structure::

    cogid	int	-	cognateid,cogid,cognateset
    entry	str	-	counterpart,word,entry,ipa,words
    taxon	str	COL	language,doculect,dialect,taxon,languages
    gloss	str	ROW	gloss,concept
    iso	        str	-	iso,isocode
    tokens	list	-	tokens,tokenized_counterpart,ipatokens

According to this structure, the first column indicates the name which is internally used to address
the given datatype. The second column indicates the program-internal datatype. The third row 
