========================================================
Sound Classes (:py:mod:`~lingpy.sequence.sound_classes`)
========================================================

.. currentmodule:: lingpy.sequence.sound_classes

This module provides functions and classes to deal with sound-class sequences. Sound classes go back
to an approach :evobib:`Dolgopolsky1964`. The basic idea behind sound classes is to reduce the IPA
space of possible phonetic segments in order to guarantee both the comparability of sounds between
languages, but also to give some assessment regarding the probability that sounds belonging to the
same class occur in correspondence relations in genetically related languages. More recently, 
sound classes have been applied in a couple of approaches, including phonetic alignment (see
:evobib:`List2012a`), and automatic cognate detection (see :evobib:`Turchin2012`,
:evobib:`List2012b`).

Functions
---------

.. autosummary::
   
   ipa2tokens
   tokens2class
   prosodic_string
   prosodic_weights
   class2tokens
   pid
   get_all_ngrams
   sampa2uni

.. currentmodule:: lingpy.data.model

Classes
-------

.. autosummary::

   Model

