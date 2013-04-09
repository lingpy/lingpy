=========
Workflows
=========

Orthographic Parsing
====================

Automatic Cognate Detection
===========================

Automatic cognate detection in LingPy is basically carried out with help of the
:py:class:`~lingpy.compare.lexstat.LexStat` class. This class takes a filename as parameter::

  >>> from lingpy import *
  >>> lex = LexStat('harry_potter.csv')

Once instantiated, a :py:class:`~lingpy.compare.lexstat.LexStat` object offers all methods and
attributes that are also defined for a :py:class:`~lingpy.basic.wordlist.Wordlist` object::

  >>> lex.width
  10
  >>> lex.height
  5



