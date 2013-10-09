Orthography Profiles
--------------------

Orthography profiles can be a combination of two files:

1. orthography profile - a document-specific specification of grapheme clusters (Unicode grapheme clusters)
2. orthography rules - a document-specific specification of orthographic replace rules specified with regular expressions

These files can be used to parse string input into orthographic tokens, as specifed in these documents. Each document should be given the proper extention and their filenames should match, e.g.:

filename.prf - for profiles
filename.rules - for rules files 

An orthography profile is often the only specification needed for orthographic tokenization. Its format is as follows:

Any number of initial lines in the file that are preceded with the hash mark "#" are used to describe the orthography profile or orthography rules files. 

# Dogon comparartive wordlist orthography profile
#
# author: Steven Moran
# author: Jelena Prokic

The first line that does not start with "#" is the header line. The header line and all following lines should contain minimally a specification of a grapheme cluster, e.g.:

graphemes
a
aa
b
ch
...

The header should always specify the keyword "graphemes". If the orthography profile author wishes to include mappings from graphemes to other elements, such as IPA representations, then columns in the orthography profile should be tab (and only tab) delimited, e.g.

# Dogon comparartive wordlist orthography profile
#
# author: Steven Moran
# author: Jelena Prokic
graphemes IPA
a	  a
aa	  aː
b	  b
ch	  tʃ
...	  ...

All columns after "graphemes" are user specified.

The other extension of the orthography profile, the orthography profile rules specification, should also contain a hash-commented header for metadata about the source document. Afterwards, orthography re-write rules are specified line-by-line as regular expressions of the follow format (NOTE: these rules are applied in order, so feeding and bleeding situations can occur):

([a|e|i|o|u])(n)(\s)([a|e|i|o|u|), \1 \2 \4

These strings are delimited by comma "," and follow the Python regular expression format. Here any vowel followed by an <n> by a space and by another vowel is replaced as first vowel followed by space followed by <n> followed by space followed by second vowel, so:

an a

becomes

a n a

Again, this is a document specific specification that should rewrite the orthographic tokenization of words, where necessary and applicable.
