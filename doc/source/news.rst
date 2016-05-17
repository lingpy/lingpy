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
:py:class:`~lingpy.align.sca.Alignments` and :py:class:`~lingpy.compare.lexstat.LexStat`. As of now, you can specify
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

For an Alignments analysis, we need to know:

* transcription,
* segments, and
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
LingPy |version| and accessible via the :py:class:`~lingpy.compare.partial.Partial` class which
itself expands the :py:class:`~lingpy.compare.lexstat.LexStat` class for automatic cognate
judgments. The main requirement for partial cognate detection is that your data is morphologically
annotated. Morphological annotation will be assumed automatically, if your dataset contains tone
annotations in South-East-Asian style, with superscript or subscript numbers, since in most of these
languages each syllable corresponds to one morpheme. If this is not given, you need to provide
morphological annotations by adding one of the currently accepted symbols for morpheme segmentation.
If you want to know which symbols are currently accepted, you can check with the
:py:class:`~lingpy.settings.rc` function::

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
  >>> part.partial_cluster(method='sca', cluster_method='upgma', ref='partial_ids')
  >>> for k in [1,2,3,4,5]:
          print(k, ''.join(part[k, 'segments']), ' '.join([str(x) for x in part[k, 'partial_ids']]))
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

  >>> alm = Alignments(test_data('partial_cognates.tsv'),ref='partial_ids', segments='segments')
  >>> alm.align()
  >>> print(SCA(alm.msa['partial_ids'][1]))
  x    u    ³⁵
  x    u    ²⁴
  x    u    ²⁴
  x    u    ²⁴

The basic idea behind this change is a different handling of etymological dictionaries within the
:py:class:`~lingpy.basic.wordlist.Wordlist` class: Earlier, an etymological dictionary was just a
transposed representation in which the cognate identifier was the key of a dictionary, and a list of values in order
of the languages linked to the words, with a zero indicating absence of values::

  >>> part = Partial(test_data('partial_cognates.tsv'), segments="segments")
  >>> etd1 = part.get_etymdict('partialids2')
  >>> etd2 = part.get_etymdict('pcogsets')
  >>> len(etd1), len(etd2)
  13, 19
  >>> etd1.keys()
  dict_keys(['7 8 9', '18 16 17', '3 4', '18 16', '13', '12', '1 1', '5 6', '1 2', '14 15', '12 10 11', '12 10', '1'])
  >>> etd2.keys()
  dict_keys([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])

In the case of ``etd1``, we reference the column "partialcogids2" for the creation of the
etymological dictionary, but for this column, there is no instruction in LingPy, no column type
defined in the wordlist-rc file. As a result, its values are interpreted as strings, and the
etymological dictionary that is created has strings as keys. In the case of the second column, the
column "pcogsets" is a reserved namespace for partial cognate sets, and all values are converted to
lists of integers. If the Wordlist class (or its descendants) creates an etymdict and detects a list
type instead of a string type or an integer type, it automatically assumes that the data is "fuzzy",
and it creates a different etymdict in which each item of the list of integer ids is transformed
into the key of the dictionary. 

Note that even if this sounds a bit confusing in the moment, we will try to further improve the
handling of this feature. This is closely related to the general change of wordlist handling in
LingPy using the wordlist-rc files, which we will try to replace with a more consistent system.

Evaluation methods
------------------

If you want to test how well a partial cognate detection analysis works, there is a special function
for evaluation, and its usage is idential with the usage of the classical
:py:class:`~lingpy.evaluate.acd.bcubes` function::

  >>> from lingpy.evaluate.acd import partial_bcubes
  >>> partial_bcubes(lex, "partialid", "pcogsets")

Furthermore, you can test n-point average precision (:py:class:`~lingpy.evaluate.acd.npoint_ap`) on string similarity measures, which we
implemented following the suggestion by :evobib:`Kondrak2002`. The n-point average precision is a
measure that can be used to test how well string similarities or distances distinguish cognates from
non-cognates, and it is therefore useful for the evaluation or the training of different string
similarity measures::

  >>> from lingpy.evaluate.acd import npoint_ap
  >>> from lingpy import *
  >>> word_pairs = [("harry", "gari"), ("walter", "woldemort"), ("harry", "walter"), ("gari", "woldemort")]
  >>> dists = [edit_dist(a, b, normalized=True) for a, b in word_pairs]
  >>> cogs = [1, 1, 0, 0] # these are our cognates for all word pairs
  >>> npoint_ap(dists, cogs)
  1.0

What's Next?
------------

.. toctree::
   :maxdepth: 1
   
   intro
   examples
   tutorial/index
   docu/index
   reference
   download
