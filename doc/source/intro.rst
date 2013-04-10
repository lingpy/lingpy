.. -*- coding: utf-8 -*-

Introduction
============


What is LingPy?
---------------

LingPy is a suite of open-source Python modules for sequence comparison,
distance analyses, data operations and visualization methods in quantitative
historical linguistics.

The main idea of LingPy is to provide a software package which, on
the one hand, integrates different methods for data analysis in quantitative
historical linguistics within a single framework, and, on the other hand,
serves as an interface for the preparation and analysis of linguistic data
using biological software packages.

What can be done with LingPy?
-----------------------------

With the help of LingPy, one can:

* tokenize and analyze IPA-encoded sequences, 
* carry out pairwise and multiple alignment analyses,
* conduct automatic cognate judgments for multiple languages, 
* calculate lexicostatistic distances between languages, 
* reconstruct language phylogenies using basic cluster algorithms, and
* export the results of these analyses to various formats which can either be
  used as input for external programs, or to visualize the results. 

How can LingPy be used?
-----------------------

LingPy is a library for the Python programming language. It can thus be used
directly from the Python shell or inside Python scripts.  The Python programming
language has many benefits. Due to its clear programming style it allows the
user to write programs very quickly while at the same time being very powerful
and easy to maintain. See http://www.python.org for details.

What is the Basic Idea Behind LingPy?
-------------------------------------

In historical linguistics, one usually deals with the following two entities
and their components:

* *signs* -- form-meaning pairs, i.e. *sequences* (words, morphemes) which
  carry a certain *functon* (meaning), and
* *languages* -- collections of *signs*.

These entities are dealt with by *comparison*. Signs and languages are compared
in different ways in order to reconstruct how they evolved from their
respective ancestor entities. Thus, when investigating the phenomena of sound
change, this is usually done by comparing the *formal* component of signs. When
comparing semantic change, the main objective is the *functional* component of
signs. And when reconstructing language phylogenies, this is done by comparing
the *signs* of different languages.

In LingPy, the different entities are represented within an object-oriented
frame work. The different aspects of linguistic signs and languages are modeled
by specific classes. Thus, if one has to deal with the formal part of the
linguistic sign, the :py:class:`~lingpy.sequence.Sequence` class in LingPy
makes it possible to deal with various automatic aspects of sound sequences,
such as, e.g., their prosodic properties, or their sound classes. For the
comparison of sign forms, one can chose between the
:py:class:`~lingpy.algorithm.classes` and the
:py:class:`~lingpy.compare.Multiple` class, which can be used to conduct
automatic alignment analyses or to calculate the similarity or distance between
sound sequences. Even language comparison is possible with help of the
:py:class:`~lingpy.lexstat.LexStat` class which can be used to conduct an
automatic search for cognates in lexicostatistical datasets or to calculate
evolutionary distances between languages.

What's Next?
-------------------------

.. toctree::
   :maxdepth: 1

   examples
   docu/index
   tutorial/index
   download

