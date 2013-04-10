=========
Workflows
=========

This is an example workflow that illustrates what can be done with LingPy. It starts from a small
dataset of the Dogon languages and leads through all important stages of a typical automatic
analysis with LingPy. 

Getting Started
===============

Make sure to have the LingPy library installed properly, see the :ref:`Installation Instructions`
for details. The dataset that will be used can be downloaded from the following link:

http://www.lingpy.org/workflow.zip

In order to start the analysis, unzip the dataset and cd into the folder. Then load the LingPy
library::

  >>> from lingpy import *

Then load the dataset into a :py:class:`~lingpy.basic.wordlist.Worlist` object by typing::
  
  >>> wl = Wordlist('DOGON.csv')

This will load the wordlist, and you can check for some specific characteristics, such as the number
of languages (the 'width' of the Wordlist), or the number of concepts (its 'height')::

  >>> wl.width
  18
  >>> wl.height
  325

Orthographic Parsing
====================

Orthographic parsing is the first step of the analysis. For the analysis, IPA format in tokenized
form is needed. In order to get the IPA tokens from the raw data, we need an orthography profile.
This has already been prepared. It is stored in the file Heath2013.prf_. In order to parse the data,
type in the following::

  >>> wl.tokenize(ortho_profile='Heath2013.prf')

As a result, LingPy will iterate over all words in the dataset and convert them to IPA, as specified
in the orthographic profile. We can now write the data to file so that we can use it in the next
step of the workflow::

  >>> wl.output('csv',filename='DOGON_tokens')

This writes the data to the file DOGON_tokens.csv_

Automatic Cognate Detection
===========================

Automatic cognate detection in LingPy is basically carried out with help of the
:py:class:`~lingpy.compare.lexstat.LexStat` class. As all classes that inherit from the :py:class`~lingpy.basic.wordlist.Wordlist` class, this class also takes a filename as parameter::

  >>> lex = LexStat('harry_potter.csv')

Once instantiated, a :py:class:`~lingpy.compare.lexstat.LexStat` object offers all methods and
attributes that are also defined for a :py:class:`~lingpy.basic.wordlist.Wordlist` object::

  >>> lex.width
  18
  >>> lex.height
  325

In order to carry out a cognate detection analysis using the LexStat method (see :evobib:`List2012b`), we first have to 
create a scoring dictionary from automatically identified sound correspondences. There are a lot of
different parameters for this analysis, but for now, the defaults should suffice, so we can call
this method by simply typing::

  >>> lex.get_scorer()

After some time that it takes to calculate all the data, we can go on with the real analysis for
cognate detection. However, since the analysis is time-consuming, it is useful to store the current
state before going on::

  >>> lex.pickle()

Now, we can carry out the cognate detection. This is a cluster method, that clusters all sequences
which are similar to each other into the same cognate set. Which sequences are clustered depends on
a threshold that we have to pass as an argument. Here we chose 0.5 as a threshold. This is a rather
conservative score that avoids that we get too many false positives::

  >>> lex.cluster(method='lexstat',threshold=0.5)

Again, we output the data. But since the LexStat method produces a lot of alternative data that is
not necessarily needed for the following analyses, we reduce the output in the CSV-format by
setting the **subset** keyword to c{True} and pass the data we want as a list to the keyword **col**. 
In order to have a nice format with all words corresponding to the same concept in the same block,
we specify the keyword **formatter** as 'concepts'::

  >>> lex.output('csv',subset=True, filename='DOGON_lexstat',subset=True,formatter='concepts',cols=['concepts','taxa','counterpart','tokens','lexstatid'])

This produces the file DOGON_lexstat.csv_ in our folder.

Phonetic Alignment
==================

Phonetic alignment is the basis of the LexStat analysis we just carried out. However, it is also
useful for the purpose of visualization, especially multiple alignment analyses can give us a quick
hint whether the cognates that an algorithm detected are "good" ones, or not. In order to carry out
multiple alignment analyses from a cognate set, we can load the data that we just wrote to file in
the previous step into an :py:class:`~lingpy.align.sca.Alignments` object. Note that we should
specify, where the cognates are. In the case of a LexStat analysis, they are stored in the
'lexstatid' column::

  >>> alm = Alignments('DOGON_lexstat.csv',cognates='lexstatid')

Carrying out a default alignment analysis is now very simple. We choose the default parameters, and
the 'library'-method for multiple alignments (see :evobib:`List2012a`), we also set the **output**
kewyord to c{True} in order to have all alignments written to separate files::

  >>> alm.align(method='library',output=True)

This will produce a new folder ``DOGON_lexstat_msa/`` that contains all multiple alignments in
separate MSA-files. More information regarding the format of these files can be found under: :ref:`msa_formats`.
The MSA-format is useful for manual editing and comparing of multiple alignments. In order to view a
whole dataset of cognate judgments and aligmnents, however, it not very appropriate. Here, LingPy
offers a specific colored HTML-output that is very helpful to inspect the results. In order to
create this output, we first have to write the aligned data to a specific format with the extension
``alm``::

  >>> alm.output('alm',cognates='lexstatid',filename='DOGON')

Now, that we created the file DOGON.alm_, we have to load the :py:func:`~lingpy.convert.plot.alm2html` from the
:py:mod:`~lingpy.convert.plot`-module. This module is not automatically loaded when importing
LingPy, so we have to import it explicitly::

  >>> from lingpy.convert.plot import alm2hmtl

Once the module is imported, we can use the function to convert the file DOGON.alm_ to colored
HTML-output::

  >>> alm2html('DOGON.alm',filename='DOGON')

As a result, we get the file DOGON.html_ in our folder.



.. _Heath2013.prf: examples/Heath2013.prf
.. _DOGON.csv: examples/DOGON.csv
.. _DOGON_tokens.csv: examples/DOGON_tokens.csv
.. _DOGON_lexstat.csv: examples/DOGON_lexstat.csv
.. _DOGON.alm: examples/DOGON.alm
.. _DOGON.html: examples/DOGON.html

