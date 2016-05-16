==================
Handling Wordlists
==================

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
all different possible **entry-types** a word can have, just as in the file `harry_potter.csv`_::

    @author: Potter, Harry
    @date: 2012-11-07
    #
    ID   CONCEPT     COUNTERPART   IPA         DOCULECT     COGID
    1    hand        Hand          hant        German       1
    2    hand        hand          hænd        English      1
    3    hand        рука          ruka        Russian      2
    4    hand        рука          ruka        Ukrainian    2
    5    leg         Bein          bain        German       3
    6    leg         leg           lɛg         English      4
    7    leg         нога          noga        Russian      5
    8    leg         нога          noha        Ukrainian    5
    9    Woldemort   Waldemar      valdemar    German       6
    10   Woldemort   Woldemort     wɔldemɔrt   English      6
    11   Woldemort   Владимир      vladimir    Russian      6
    12   Woldemort   Володимир     volodimir   Ukrainian    6
    13    Harry       Harald        haralt      German       7
    14   Harry       Harry         hæri        English      7
    15   Harry       Гарри         gari        Russian      7
    16   Harry       Гаррi         hari        Ukrainian    7

This format is, of course, much more redundant, than the word list format, but
it allows to display multiple entry-types for the counterparts of a given
concept in a given language. Moreover, this format is the basic of the
:py:class:`~lingpy.basic.wordlist.Wordlist` class in LingPy, which makes it easy
to handle word lists with multiple entry-types of words.

Basic Operations with Help of Wordlists
---------------------------------------

The above-given csv-file `harry_potter.csv`_ is available in the test folder of LingPy.
In order to get it loaded, we simply pass it as first argument to the Wordlist
class::
    
    >>> from lingpy import *
    >>> d = Wordlist('harry_potter.csv')

If one wants to access only the IPA values in tabular format, all one has to do
is::

    >>> wl.ipa
    [['wɔldemɔrt', 'valdemar', 'vladimir', 'volodimir'],
     ['hæri', 'haralt', 'gari', 'hari'],
     ['lɛg', 'bain', 'noga', 'noha'],
     ['hænd', 'hant', 'ruka', 'ruka']]

The same for cognates::

    >>> wl.cognate
    [[6, 6, 6, 6], [7, 7, 7, 7], [4, 3, 5, 5], [1, 1, 2, 2]]

Or for the languages and the concepts in the dataset::

    >>> wl.language
    ['English', 'German', 'Russian', 'Ukrainian']
    >>> wl.concept
    ['Harry', 'Woldemort', 'hand', 'leg']
    
Furthermore, using specific functions, even more concise samples of the data can
be extracted, thus, using the
:py:class:`~lingpy.basic.wordlist.Wordlist.get_dict` function, we can specify a
given language and extract all phonetic transcriptions corresponding to a given
concept as a dictionary::

    >>> wl.get_dict(col="German",entry="ipa")
    {'Harry': ['haralt'],
     'Woldemort': ['valdemar'],
     'hand': ['hant'],
     'leg': ['bain']}

We can likewise extract all cognate IDs corresponding to a given concept by
using the function :py:class:`~lingpy.basic.wordlist.Wordlist.get_list`::

    >>> wl.get_list(row="hand",entry="cogid",flat=True)
    [1, 1, 2, 2]
    
Other entry-types can be added::

    >>> from lingpy.algorithm.misc import ipa2tokens
    >>> wl.add_entries("tokens","ipa",ipa2tokens)
    >>> wl.tokens
    [[['w', 'ɔ', 'l', 'd', 'e', 'm', 'ɔ', 'r', 't'],
      ['v', 'a', 'l', 'd', 'e', 'm', 'a', 'r'],
      ['v', 'l', 'a', 'd', 'i', 'm', 'i', 'r'],
      ['v', 'o', 'l', 'o', 'd', 'i', 'm', 'i', 'r']],
     [['l', 'ɛ', 'g'],
      ['b', 'ai', 'n'],
      ['n', 'o', 'g', 'a'],
      ['n', 'o', 'h', 'a']],
     [['h', 'æ', 'n', 'd'],
      ['h', 'a', 'n', 't'],
      ['r', 'u', 'k', 'a'],
      ['r', 'u', 'k', 'a']],
     [['h', 'æ', 'r', 'i'],
      ['h', 'a', 'r', 'a', 'l', 't'],
      ['g', 'a', 'r', 'i'],
      ['h', 'a', 'r', 'i']]]
    
The wordlist.rc file
--------------------

The structure of word lists is defined by the configuration file `wordlist.rc`_. This file is
automatically loaded when initializing a Wordlist instance::

    >>> wl = Wordlist(data)

It can, however, also be passed by the user::

    >>> wl = Wordlist(data, conf="path_to_file")

All rc-files (which are used for different wordlist-like object in LingPy)
are currently located at `lingpy/data/rc/` and have
a simple tab-separated structure of four three columns:

1. basic namespace (alphanumeric, lower case)
2. class of all entries in that namespace
3. alias for the namespace (alphanumeric, all lower case, comma-separated)

As an example, consider the following minimal rc-file for a
wordlist object::

   ipa     str    ipa,orthography,transcription
   tokens  lambda x: x.split(" ")  tokens,segments
   cogid   int     cognate_set_id,cognates,cogid

This rc-file, which you can call by passing the path of your file as an
argument when loading a wordlist, will treat all entries in columns named
"ipa" or "orthography" or "transcription" in your data as strings, it will
further define "ipa" as the basic name for those columns and use this name when
you output the file. It will split all entries in the column "tokens" (or
"segments") along spaces and store them as a list, and it will convert all
"cogid" entries to integers.

.. _harry_potter.csv: examples/harry_potter.csv
.. _wordlist.rc: examples/wordlist.rc
