Helper Functions for SCA Alignment (:py:class:`~lingpy.algorithm.cython.calign`, and :py:class:`~lingpy.algorith.cython.misc`)
==============================================================================================================================

The helper functions and classes below play an important role in all SCA alignment algorithms in LingPy (:evobib:`List2012b`).
They are implemented both in pure Python and in Cython (only supported for Python 3), in order to
allow for faster implementations of the core alignment functions. Instead of using these functions
directly, we recommend to use the more general functions which you can find in the
:py:class:`~lingpy.align.pairwise` and the :py:class:`~lingpy.align.multiple` module of LingPy, and
which are based on the helper functions we list below.

.. currentmodule:: lingpy.algorithm.cython.calign

Functions
---------

.. autosummary::

   ~lingpy.algorithm.cython.calign.globalign
   ~lingpy.algorithm.cython.calign.secondary_globalign
   ~lingpy.algorithm.cython.calign.localign
   ~lingpy.algorithm.cython.calign.secondary_localign
   ~lingpy.algorithm.cython.calign.semi_globalign
   ~lingpy.algorithm.cython.calign.secondary_semi_globalign
   ~lingpy.algorithm.cython.calign.dialign
   ~lingpy.algorithm.cython.calign.secondary_dialign
   ~lingpy.algorithm.cython.calign.align_pair
   ~lingpy.algorithm.cython.calign.align_pairwise
   ~lingpy.algorithm.cython.calign.align_pairs
   ~lingpy.algorithm.cython.calign.align_profile
   ~lingpy.algorithm.cython.calign.score_profile
   ~lingpy.algorithm.cython.calign.swap_score_profile
   ~lingpy.algorithm.cython.calign.corrdist

Classes
-------

.. autosummary::

   ~lingpy.algorithm.cython.misc.ScoreDict
  

