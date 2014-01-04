=============
Basic Formats
=============

In the following, we list some of the formats that are frequently used by LingPy, be it that they
are taken as input formats, or that they are produced as output from the classes and methods
provided by LingPy.

.. _wordlist_format:

Wordlist-Format: Basic Format for Storing Large Datasets
--------------------------------------------------------

For the :py:class:`~lingpy.basic.wordlist.Wordlist` class (and also for all classes that inherit
from it, such as :py:class:`~lingpy.compare.lexstat.LexStat`,
:py:class:`~lingpy.compare.phylogney.PhyBo`, :py:class:`~lingpy.align.sca.Alignments`), a
simple csv-format is used. This format is a simple tab-delimited text file in which the header
specifies all entry types in a given dataset::

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
    9    Harry       Harald        haralt      German       7
    10   Harry       Harry         hæri        English      7
    11   Harry       Гарри         gari        Russian      7
    12   Harry       Гаррi         hari        Ukrainian    7

This format can be further extended by adding key-value pairs in the lines before the header, such
as, for example, information regarding the author, the data, or general notes:: 

    @author: Potter, Harry
    @date: 2012-11-07
    @note: Be careful with this data, it might have been charmed...
    # 
    ID   CONCEPT     COUNTERPART   IPA         DOCULECT     COGID
    1    hand        Hand          hant        German       1
    2    hand        hand          hænd        English      1
    3    hand        рука          ruka        Russian      2
    ...  ...         ...           ...         ...          ...

This format is, of course, rather redundant, but it allows to display multiple entry-types for
language data. Furthermore, the data can be easily extended. Thus, one can add multiple alignments,
using the standard formats for multiple alignments, as described under :ref:`msa_formats`,
by enclosing them in specific html-tags and placing them before the real data::
  
    @author: Potter, Harry
    @date: 2012-11-07
    @note: Be careful with this data, it might have been charmed...
    # 
    <msa id="6" ref="cogid">
         Harry Potter Testset
         Woldemort (in different languages)
         English     w    o    l    -    d    e    m    o    r    t
         German.     w    a    l    -    d    e    m    a    r    -
         Russian     v    -    l    a    d    i    m    i    r    -
         Ukrainian   v    o    l    o    d    y    m    y    r    -
    </msa>
    # 
    ID   CONCEPT     COUNTERPART   IPA         DOCULECT     COGID
    1    hand        Hand          hant        German       1
    2    hand        hand          hænd        English      1
    3    hand        рука          ruka        Russian      2
    ...  ...         ...           ...         ...          ...

.. _alignment_formats:

Basic Formats for Phonetic Alignments
-------------------------------------

.. _psa_formats:

Pairwise Alignments (PSQ and PSA)
---------------------------------

The input format for text files containing unaligned sequence pairs is called PSQ-format. Files
in this format should have the extension ``psq``. The first line of a PSQ-file contains information
regarding the dataset. The sequence pairs are given in triplets, with a sequence identifier in the
first line of a triplet (containing the meaning, or orthographical information) and the two
sequences in the second and third line, whereas the first column of each sequence line contains the
name of the taxon and the second column the sequence in IPA format. All triplets are divided by one
empty line. As an example, consider the file `harry_potter.psq`_: 

.. raw:: html
  :file: examples/harry_potter.psq.html

The output counterpart of the PSQ-format is the PSA-format. It is a specific format for text files
containing already aligned sequence pairs. Files in this format should have the extension ``psa``. The
first line of a PSA-file contains information regarding the dataset. The sequence pairs are given in
triplets, with a sequence identifier in the first line of a triplet (containing the meaning, or
orthographical information) and the aligned sequences in the second and third line, whith the name
of the taxon in the first column and all aligned segments in the following columns, separated by
tabstops. All triplets are divided by one empty line. As an example, consider the file
harry_potter.psa_:

.. raw:: html
  :file: examples/harry_potter.psa.html

.. _msa_formats:

Multiple Alignments (MSQ and MSA)
---------------------------------

A specific format for text files containing multiple unaligned sequences is the MSQ-format.
Files in this
format should have the extension ``msq``. The first line of an msq-file contains information regarding
the dataset. The second line contains information regarding the sequence (meaning, identifier), and
the following lines contain the name of the taxa in the first column and the sequences in IPA format
in the second column, separated by a tabstop. As an example, consider the file harry_potter.msq_:

.. raw:: html
  :file: examples/harry_potter.msq.html

The msa-format is a specific format for text files containing already aligned sequence pairs. Files
in this format should have the extension ``msa``. The first line of a MSA-file contains information
regarding the dataset. The second line contains information regarding the sequence (its meaning, the
protoform corresponding to the cognate set, etc.). The aligned sequences are given in the following
lines, whereas the taxa are given in the first column and the aligned segments in the following
columns. Additionally, there may be a specific line indicating the presence of swaps and a specific
line indicating highly consistent sites (local peaks) in the MSA. The line for swaps starts with the
headword SWAPS whereas a plus character (+) marks the beginning of a swapped region, the dash
character (-) its center and another plus character the end. All sites which are not affected by
swaps contain a dot. The line for local peaks starts with the headword LOCAL. All sites which are
highly consistent are marked with an asterisk (*), all other sites are marked with a dot (.). As an
example, consider the file harry_potter.msa_:

.. raw:: html
  :file: examples/harry_potter.msa.html


.. _csv_alignment_format:

.. _harry_potter.psq: examples/harry_potter.psq
.. _harry_potter.psa: examples/harry_potter.psa
.. _harry_potter.msq: examples/harry_potter.msq
.. _harry_potter.msa: examples/harry_potter.msa
