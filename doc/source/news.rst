===============================
Whats New in Version |version|?
===============================

Version |version| of LingPy introduces a couple of new algorithms, but more importantly, much more
consistency than earlier versions. As of now, LingPy runs in Python 2 and Python 3 across all
major platforms, and we try hard to make sure it stays in this way. 

Below, we list some interesting new features shipped along with LingPy |version| and which you might
want to test.

Enhanced Documentation
----------------------

You may not notice this directly, but we have enhanced the documentation. New algorithms which have
been added are all fully documented, but in addition, we have finally updated the documentation for
the deep alignment modules, such as :py:class:`~lingpy.algorithm.cython.calign` and
:py:class:`~lingpy.algorith.cython.misc`. If you want to use these functions to speed up your
alignments, or to create your own algorithms that build on LingPy's deep alignment functions, you
should definitely have a look at the additional documentation and some examples which we added.

The Command Line Interface
--------------------------

Version |version| introduces a command line interface which remains purely experimental for the
moment, but will be further enhanced in the future in order to enable users to include LingPy in
their pipelines with other programs.

As an example, after installing LingPy on your machine, just open a terminal, and try the
following::

  $ lingpy -h

This will show you all basic instructions regarding the usage of the command line interface. The
basic idea here is that we split the commandline into different subtasks::
  
  $ lingpy pairwise -h
  usage: lingpy pairwise [-h] [--input-file INPUT_FILE]
                       [--output-file OUTPUT_FILE] [--factor FACTOR]
                       [--gop GOP] [--scale SCALE]
                       [--restricted-chars RESTRICTED_CHARS]
                       [--strings STRINGS STRINGS]
                       [--mode {global,local,overlap,dialign}] [--distance]
                       [--method {sca,basic}]
  ...

If you want to align a couple of strings, you can, for example, do the following::

  $ lingpy multiple -s woldemort walter waldemar
  w    o    l    d    e    m    o    r    t
  w    a    l    t    e    -    -    r    -
  w    a    l    d    e    m    a    r    -

There are many more possibilities, but since this features is still experimental, we did not work
out a full documentation, so for the moment, we recommend to follow the instructions provided by
adding the "-h" command to each of the subcommands currently available, when testing this feature.

Slowly Saying Goodbye to the Alias System
-----------------------------------------

The alias-system in LingPy wordlists has lead to some confusion among users. Basically, the idea was
to guarantee a maximum of flexibility. As a result, the flexibility lead to a high degree of
confusion. We should note, however, that most of these aspects are documented, especially in our
tutorial section, and it is recommended for all users and those interested in using LingPy to give
these a proper read.

As of version |version|, however, we re-arranged the handling to potentially overcome the problem of
namespaces by adding explicit arguments to the main classes for alignments and cognate sets, like
:py:class:`~lingpy.align.sca.Alignments` and :py:class:`~lingpy.compare.lexstat.LexStat.` As of now, you can specify
integral parts of the data by passing how you define your namespace as an argument upon
initialization. For a LexStat analysis, for example, we need to know where we find the following
information and how it should be named:

* transcription (usually called "ipa")
* segments (usually called "tokens", if not given, it is created from "transcription")
* classes (the sound classes, which are usually created from "segments")
* langid (language identifiers, which are needed to create the individual segments for each language, usually created)
* prostrings (the prosodic strings, distinguishing context, usually created)
* numbers (the combination of langid, classes, and prostrings to form unique segment representations for each segment in each language, also the basic value passed to the scoring function)
* duplicates (the column in which duplicates are stored, and marked by a 1 for a duplicated entry, that is, an entry which is identical with another entry with a different meaning, usually created)
* alignment (for alignment analyses, the argument where alignments should be stored, created automatically or manually)

From now on, you can define your own namespaces by passing those as arguments when loading a
wordlist for a LexStat or an alignment analysis::

   >>> lex = LexStat('myfile.tsv', segments="segments", transcription="transcription")
   >>> alm = Alignments('myalms.tsv', segments="segments", alignment="myperfectalignment")

This allows for a much greated flexibility and is a first step towards replacing the alias and the
rc-configuration file system in which the namespace was handled previously. 

Partial Cognate Detection
-------------------------

We present an initial attempt to carry out partial cognate detection and which is shipped along with
LingPy |version| and accessible via the :py:class:`~lingpy.compare.partial.Partial` class which itself expands
the :py:class:`~lingpy.compare.lexstat.LexStat` class for automatic cognate judgments. The main requirement for
partial cognate detection is that your data is morphologically annotated. Morphological annotation
will be assumed automatically, if your dataset contains tone annotations in South-East-Asian style,
with superscript or subscript numbers, since in most of these languages each syllable corresponds to
one morpheme. If this is not given, you need to provide morphological annotations by adding one of
the currently accepted symbols for morpheme segmentation. If you want to know which symbols are
currently accepted, you can check with the :py:class:`~lingpy.settings.rc` function::

  >>> from lingpy import rc
  >>> rc('morpheme_separators')
  '◦+'
  >>> rc('word_separators')
  '_#'

As you can see, we currently support two morpheme separators and two word separators. Note, however,
that there is no difference in treatment when using the new method for partial cognate detection:
whether your data uses the morpheme separators or the word separators should not make a difference,
since LingPy will split all words in your data into their smallest units.

In order to test partial cognate detection, you can use a very small test set shipped with the
LingPy test suite::

  >>> from lingpy.compare.partial import Partial
  >>> from lingpy.tests.util import test_data
  >>> part = Partial(test_data('partial_cognates.tsv'), segments="segments")
  >>> part.partial_cluster(method='sca', cluster_method='upgma', ref='pcogs')
  >>> for k in [1,2,3,4,5]:
          print(k, ''.join(part[k, 'segments']), ' '.join([str(x) for x in part[k, 'pcogs']]))
  ...
  1 xu²⁴+ni⁵⁵ 3 4
  2 xu²⁴ni⁴⁴ 3 4
  3 xu³⁵+ni⁵⁵ 3 4
  4 bu¹³li⁵³ 1 2
  5 bu¹³ 1

Note that now, if you provide both partial cognate sets and segmented words, you can even align all
partial cognate sets in separation. The only current restriction here is that you need to follow our
namespace for fuzzy cognates as defined in the `wordlist.rc` (but you can modify it accordingly by
just creating your own `wordlist.rc` and passing it as an extra argument)::

  >>> alm = Alignments(test_data('partial_cognates.tsv'),ref='partial_cognate_sets', segments='segments')
  >>> alm.align()
  >>> print(SCA(alm.msa['pcogsets'][1]))
  x    u    ³⁵
  x    u    ²⁴
  x    u    ²⁴
  x    u    ²⁴


What's Next?
------------

.. toctree::
   :maxdepth: 1
   
   intro
   examples
   tutorial/index
   docu/index
   reference/modules
   download
