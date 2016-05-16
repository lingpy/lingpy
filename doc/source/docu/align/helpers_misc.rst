Helper Functions for Traditional Alignment (:py:class:`~lingpy.algorithm.cython.talign`)
========================================================================================

The helper functions and classes below play an important role in traditional alignment algorithms in
LingPy which do not make use of sound classes.  They are implemented both in pure Python and in
Cython (only supported for Python 3), in order to allow for faster implementations of the core
alignment functions. Instead of using these functions directly, we recommend to use the more general
functions which you can find in the :py:class:`~lingpy.align.pairwise` and the
:py:class:`~lingpy.align.multiple` module of LingPy, and which are based on the helper functions we
list below.

.. currentmodule:: lingpy.algorithm.cython.calign

Functions
---------

.. autosummary::

   ~lingpy.algorithm.cython.talign.globalign
   ~lingpy.algorithm.cython.talign.localign
   ~lingpy.algorithm.cython.talign.semi_globalign
   ~lingpy.algorithm.cython.talign.dialign
   ~lingpy.algorithm.cython.talign.align_pair
   ~lingpy.algorithm.cython.talign.align_pairwise
   ~lingpy.algorithm.cython.talign.align_pairs
   ~lingpy.algorithm.cython.talign.align_profile
   ~lingpy.algorithm.cython.talign.score_profile
   ~lingpy.algorithm.cython.talign.swap_score_profile

  

