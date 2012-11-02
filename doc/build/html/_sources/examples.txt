========
Examples
========

The Very Basics
---------------

Load the library::

    >>> import lingpy

Retrieve basic information about the library::

    >>> print(lingpy.__doc__)
    
    LingPy --- A Python library for quantitative historical linguistics
    ===================================================================
    
    Documentation is available in the docstrings. Online documentation is available
    at http://lingulist.de/lingpy/
    
    Submodules
    -----------
    sequence                --- Sequence Modelling
    compare                 --- Sequence Comparison
    lexstat                 --- Language Comparison
    
    Subpackages
    -----------
    algorithm               --- Basic Algorithms for Sequence Comparison
    align                   --- Specific Algorithms for PSA and MSA
    data                    --- Data Handling
    output                  --- Output Handling
    test                    --- Tests and Evaluation

Load all important packages::

    >>> from lingpy import *

Pairwise Alignment Analyses
---------------------------

Load the data given in the file `test.psq <http://lingulist.de/lingpy/test.psq.txt>`_::

    >>> pairs = Pairwise(get_file('test.psq'))

Check for the sequence pairs::

    >>> pairs.seqs
    [('waldemar', 'vladimir'), ('woldemort', 'vladimir'), ('woldemort', 'waldemar')]

Align all sequences pairwise::

    >>> pairs.align()

Print the data to the screen::

    >>> print(pairs)
    w    a    l    -    d    e    m    a    r
    v    -    l    a    d    i    m    i    r
    31.0
    
    w    a    l    -    d    e    m    a    r
    v    -    l    a    d    i    m    i    r
    31.0
    
    w    o    l    -    d    e    m    o    r    t
    v    -    l    a    d    i    m    i    -    r
    31.0
    
    w    o    l    d    e    m    o    r    t
    w    a    l    d    e    m    a    -    r
    53.0


Multiple Alignment Analyses
---------------------------

Load the file `test.msq <http://lingulist.de/lingpy/test.msq.txt>`_::

    >>> mult = Multiple(get_file('test.msq'))

Align the data using a simple progressive algorithm::
    >>> mult.prog_align()

Print the data to the screen::

    >>> print mult
    w    -    o    l    d    e    m    o    r    t
    w    -    a    l    d    e    m    a    -    r
    v    l    a    -    d    i    m    i    -    r

Align the data using the library method::

    >>> mult.lib_align()
    >>> print mult
    w    o    l    -    d    e    m    o    r    t
    w    a    l    -    d    e    m    a    -    r
    v    -    l    a    d    i    m    i    -    r

Carry out a check for swapped sites::

    >>> mult.swap_check()
    True

Get the percentage identity of the data::

    >>> mult.get_pid()
    0.43333333333333335

Get the sum-of-pairs score of the data::

    >>> mult.sum_of_pairs()
    7.2166666666666668

Automatic Search for Cognates
-----------------------------

Load the file `SLV.lxs <http://lingulist.de/lingpy/SLV.lxs.txt>`_ from the test sets::

    >>> lex = LexStat(get_file('SLV.lxs'))

Conduct an automatic search for cognates::

    >>> lex.analyze(0.6)
    [i] Loaded and calculated all essential values.
    [i] Calculating scores for Russian and Russian ...
    [i] Calculating scores for Russian and Polish ...
    [i] Calculating scores for Russian and Bulgarian ...
    [i] Calculating scores for Russian and Czech ...
    [i] Calculating scores for Polish and Polish ...
    [i] Calculating scores for Polish and Bulgarian ...
    [i] Calculating scores for Polish and Czech ...
    [i] Calculating scores for Bulgarian and Bulgarian ...
    [i] Calculating scores for Bulgarian and Czech ...
    [i] Calculating scores for Czech and Czech ...
    [i] Created the library.
    [i] Calculated pairwise scores.
    [i] Calculated cognates.

Write the results of the analysis to the file `SLV.alm <http://lingulist.de/lingpy/SLV.alm.txt>`_ (aligned
output format)::

    >>> lex.output('alm')

Convert the text-based ``alm``-format into colored HTML (see `SLV.html <http://lingulist.de/lingpy/SLV.html>`_)::

    >>> alm2html('SLV.alm')

Calculate the pairwise distances::

    >>> lex.pairwise_distances()
    array([[ 0.        ,  0.21818182,  0.28181818,  0.24770642],
           [ 0.21818182,  0.        ,  0.29090909,  0.1559633 ],
           [ 0.28181818,  0.29090909,  0.        ,  0.31192661],
           [ 0.24770642,  0.1559633 ,  0.31192661,  0.        ]])

Cluster the data using the Neighbor-Joining algorithm::

    >>> neighbor(lex.pairwise_distances(),lex.taxa)
    '(((Russian:0.11,Bulgarian:0.18):0.05,Polish:0.07):0.09,Czech:0.09);'

Cluster the data using the UPGMA algorithm::

    >>> upgma(lex.pairwise_distances(),lex.taxa)
    '(Bulgarian:0.15,(Russian:0.12,(Polish:0.08,Czech:0.08):0.12):0.15);'

What's Next?
------------

.. toctree::
   :maxdepth: 1

   docu/index
   download
