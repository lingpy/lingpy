===============================================
Word Lists (:py:class:`~lingpy.basic.wordlist`)
===============================================

Word lists represent the core of LingPy's data model, and a proper understanding of how we deal with
word lists is important for automatic cognate detection, alignments, and borrowing detection. The
basic class that handles word lists, is the :py:class:`~lingpy.basic.wordlist.Wordlist` class, which
is also the base class of the :py:class:`~lingpy.compare.lexstat.LexStat` class for automatic
cognate detection and the :py:class:`~lingpy.align.sca.Alignments` class for multiple alignment of
cognate words.

.. currentmodule:: lingpy.basic.wordlist

Functions
---------

.. autosummary::

   ~lingpy.basic.wordlist.get_wordlist

Classes
-------

.. autosummary:: 

   ~lingpy.basic.wordlist.Wordlist

