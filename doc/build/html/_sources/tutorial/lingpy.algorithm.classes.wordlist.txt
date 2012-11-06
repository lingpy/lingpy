==================
The WordList Class
==================

What is a Word List?
--------------------

Generally, a word list is a simple tabular data structure in which multiple
languages are structured in such a way that words are ordered in rows and
columns according to the language to which they belong and the concept they
denote. 

How are Word Lists defined?
---------------------------

In LingPy the WordList class handles wordlists. 

The wordlist.conf file
----------------------

The structure of word lists is defined by the configuration file `wordlists.conf`. This file is
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
