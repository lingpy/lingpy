========
Examples
========

The Very Basics
---------------

Load the library::

    >>> import lingpy

    Retrieve basic information about the library::
    
    >>> print(lingpy.__doc__)
    LingPy package for quantitative tasks in historical linguistics.
    
    Documentation is available in the docstrings. Online documentation is available
    at http://lingpy.org
    
    Subpackages
    -----------
    algorithm  --- Basic Algorithms for Sequence Comparison
    align      --- Specific Algorithms Alignment Analyses
    basic      --- Basic Classes for Language Comparison
    check      --- Classes for Exceptions, Warnings, and Check
    compare    --- Basic Modules for Language Comparison
    convert    --- Functions for Format Conversion
    data       --- Data Handling
    evaluate   --- Basic Classes and Functions for Algorithm Evaluation
    read       --- Basic Functions for Data Input
    sequence   --- Basic Functions for Sequence Modeling
    thirdparty --- Temporary Forks of Third-Party-Modules


Load all important packages::

    >>> from lingpy import *

Alignment Analyses
------------------

Carry out a multiple alignment analysis of four sequences. First, define the sequences::

   >>> seqs = ['woldemort','waldemar','wladimir','vladymyr']

Create an instance of the :py:class:`~lingpy.align.multiple.Multiple` class::

   >>> msa = Multiple(seqs)

Align the sequences, using simple progressive alignment with default parameters::

   >>> msa.prog_align()

Print out the results to the screen::

   >>> print(msa)
   w    o    l    -    d    e    m    o    r    t
   w    a    l    -    d    e    m    a    r    -
   w    -    l    a    d    i    m    i    r    -
   v    -    l    a    d    y    m    y    r    -

Reconstruction of Phylogenetic Trees
------------------------------------

Import LingPy's cluster module::

  >>> from lingpy.algorithm import cluster

Set up a couple of languages::

  >>> languages = ['Norwegian','Swedish','Icelandic','Dutch','English']
  
Define a distance matrix::

  >>> distances = cluster.squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])

Carry out a Neighbor-Joining cluster analysis::

  >>> tree = neighbor(distances,languages)
  >>> print(tree)
  '(((Norwegian:0.18,(Swedish:0.12,Icelandic:0.28):0.21):0.17,Dutch:0.31):-0.01,English:-0.01);'

Print the tree in ASCII-Art to the screen (using PyCogent's tree class)::

  >>> tree = LoadTree(treestring=tree)
  >>> print(tree.asciiArt())
                                /-Norwegian
                      /edge.1--|
                     |         |          /-Swedish
            /edge.2--|          \edge.0--|
           |         |                    \-Icelandic
  -root----|         |
           |          \-Dutch
           |
            \-English
  
  
What's Next?
------------

.. toctree::
   :maxdepth: 1
   
   tutorial/index
   docu/index
   download
