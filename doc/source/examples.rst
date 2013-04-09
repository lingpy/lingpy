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

Align the sequences, using simpole progressive alignment with default parameters::

   >>> msa.prog_align()

Print out the results to the screen::

   >>> print(msa)
   w    o    l    -    d    e    m    o    r    t
   w    a    l    -    d    e    m    a    r    -
   w    -    l    a    d    i    m    i    r    -
   v    -    l    a    d    y    m    y    r    -


What's Next?
------------

.. toctree::
   :maxdepth: 1

   docu/index
   tutorial/index
   download
