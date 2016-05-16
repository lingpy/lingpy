Miscellaneous Helper Functions (:py:class:`~lingpy.algorithm.cython.malign`)
============================================================================

The helper functions below are miscellaneous "deep" implementations of alignment and string
similarity algorithms. 
They are implemented both in pure Python and in
Cython (only supported for Python 3), in order to allow for faster implementations of the core
alignment functions. Instead of using these functions directly, we recommend to use the more general
functions which you can find in the :py:class:`~lingpy.align.pairwise` and the
:py:class:`~lingpy.align.multiple` module of LingPy, and which are based on the helper functions we
list below.

.. currentmodule:: lingpy.algorithm.cython.calign

Functions
---------

.. autosummary::

   ~lingpy.algorithm.cython.malign.edit_dist
   ~lingpy.algorithm.cython.malign.nw_align
   ~lingpy.algorithm.cython.malign.restricted_edit_dist
   ~lingpy.algorithm.cython.malign.structalign
   ~lingpy.algorithm.cython.malign.sw_align
   ~lingpy.algorithm.cython.malign.we_align

  

