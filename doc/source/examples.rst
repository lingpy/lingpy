========
Examples
========

The Very Basics
---------------

Load the library::

    >>> import lingpy
    >>> import lingpy as lp

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

Use a short-cut to align the four sequences:

   >>> mult_align(['woldemort', 'waldemar', 'wladimir', 'vladymyr'], pprint=True)

Reconstruction of Phylogenetic Trees
------------------------------------

Import LingPy::

  >>> from lingpy import *

Import the squareform-function for convenience::

  >>> from lingpy.algorithm import squareform

Set up a couple of languages::

  >>> languages = ['Norwegian','Swedish','Icelandic','Dutch','English']
  
Define a distance matrix::

  >>> distances = squareform([0.5,0.67,0.8,0.2,0.4,0.7,0.6,0.8,0.8,0.3])

Carry out a Neighbor-Joining cluster analysis::

  >>> tree = neighbor(distances,languages)
  >>> print(tree)
  '(((Norwegian:0.18,(Swedish:0.12,Icelandic:0.28):0.21):0.17,Dutch:0.31):-0.01,English:-0.01);'

Print the tree in ASCII-Art to the screen (using PyCogent's tree class)::

  >>> tree = Tree(tree)
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
   reference
   download
