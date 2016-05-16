================
Workflow Example
================

This is an example workflow that illustrates some of the functionality of LingPy. We will use a
small dataset by :evobib:`Kessler2001` in order to illustrate how to get from word list data to
aligned cognate sets.

We start by loading the data, which is located in LingPy's test suite, and can be accessed with help
of the :py:class:`~lingpy.tests.util.test_data` function. We use this function to load the file `KSL.qlc` (it provides the exact path to the file), and load the file as a :py:class:`~lingpy.compare.lexstat.LexStat` object::

  >>> from lingpy import *
  >>> from lingpy.tests.util import test_data
  >>> lex = LexStat(test_data('KSL.qlc')

After we loaded the data, which is given in the general LingPy format for wordlists (see
:doc:`lingpy.basic.wordlist` for details), we can access some of its basic statistics, which are
provided by the :py:class:`~lingpy.basic.wordlist.Wordlist` class upon which the LexStat class is
based, like width (number of languages), or height (number of concepts), or length (number of words)::

  >>> lex.width, lex.height, len(lex)
  (7, 200, 1400)
  >>> lex.cols
  ['Albanian', 'English', 'French', 'German', 'Hawaiian', 'Navajo', 'Turkish']

We can also determine the coverage, which is a dictionary with language name as key and number of
concepts as value::

  >>> lex.coverage()
  {'Albanian': 200,
   'English': 200,
   'French': 200,
   'German': 200,
   'Hawaiian': 200,
   'Navajo': 200,
   'Turkish': 200}

Let's start and search for cognates now. If all works fine, and your data is in plain IPA (without
any strange characters and erros), this should work out of the box. However, if you encounter
difficulties, we recommend to make an explicit check by loading the LexStat object with the
check-keyword set to **True**::

  >>> lex = LexStat(test_data('KSL.qlc', check=True)

If the class is loaded normally, this means, your data is fine. But now let's start with cognate
judgments, and let's try the LexStat algorithm (:evobib:`List2014d`), which derives individual scorers for all language
pairs based on randomized alignments. In order to use LexStat, we first need to compute the scorer. There are many parameters, but let's just stick to defaults for now::  
  
  >>> lex.get_scorer()

Once this is done, we can compute the cognate sets. Here, we should tell LingPy in which column they
should be stored, so we pass the keyword argument "ref" and specify it as "cognates". We also need
to pass a threshold. Here, for the LexStat method, a threshold of 0.6 normally works quite well::

  >>> lex.cluster(method="lexstat", threshold=0.6, ref="cognates")

Now, we could already write the data to file, where we choose "tsv" as fileformat and specify
"KSL_new" as our filename::

  >>> lex.output('tsv', filename="KSL_new")

This resulting output file will be very large, since it contains all parameters and the computed
scoring functions for all language pairs. If you want to avoid that and only have the raw TSV
format, specify the "ignore" keyword and set "prettify" to "False"::

  >>> lex.output('tsv', filename="KSL_new", ignore="all", prettify=False)

How well was this automatic cognate detection? Let's test it by computing the B-Cubed scores
(:evobib:`Bagga1999`). These rank between 0 and 1, with one being good and 0 being bad, and come in
three flavors of precision (amount of false positives), recall (amount of false negatives) and
F-Score (combined score). We compute them by passing the wordlist object, the gold standard (stored
in column "cogid" in our "KSL.qlc" file), and our computed cognate sets (stored in "cognates", as we
specified using the "ref" keyword)::

  >>> from lingpy.evaluate.acd import bcubes, diff
  >>> bcubes(lex, "cogid", "cognates")
  *************************'
  * B-Cubed-Scores        *
  * --------------------- *
  * Precision:     0.9293 *
  * Recall:        0.9436 *
  * F-Scores:      0.9364 *
  *************************'
  (0.9292857142857142, 0.9435714285714284, 0.9363740873923939)

As we can see, the score is rather high, but we should keep in mind that the dataset is rather
small. Note that the scores may differ on your computer, since LexStat shuffles the word lists to
create a random distribution. Normally, however, the differences should not be too huge.

Now that we know that our data has been properly analyzed with good cognate scores, lets align it,
using the :py:class:`~lingpy.align.sca.Alignments` class. We can directly initialize it from the
LexStat object, but we need to pass the "cognates" column as "ref" (this tells LingPy, where to
search for cognate sets which are then multiply aligned), and then, we call the function
:py:class:`~lingpy.align.sca.Alignments.align`, using the defaults for convenience::

  >>> alm = Alignments(lex, ref='cognates')
  >>> alm.align()

If you want to see the results of this analysis, you need to write them to file. But there, it is
also difficult to see the alignments, since they are in a TSV-file in just another column, called
"alignment" as a default. Another possibility is to write data to HTML format instead. This means
you can't re-import the data into LingPy, but you can inspect the results::

  >>> alm.output('html', filename="KSL")

This will create the file KSL.html_ which you can inspect by loading it in your webbrowser.

Finally, let's make a tree of the data. This is very straightforward by passing the "cognates"
column as a reference to the :py:class:`~lingpy.basic.wordlist.Wordlist.calculate` function. Printing of the resulting "tree" which is
created as an attribute of the LexStat object, is possible with help of LingPy's
:py:class:`~lingpy.basic.tree.Tree` class::

  >>> lex.calculate('tree', ref='cognates')
  >>> print(lex.tree.asciiArt())
            /-Hawaiian
           |
  -root----|          /-Turkish
           |         |
            \edge.4--|          /-Navajo
                     |         |
                      \edge.3--|                    /-English
                               |          /edge.0--|
                               |         |          \-German
                                \edge.2--|
                                         |          /-Albanian
                                          \edge.1--|
                                                    \-French
  
Well, given that there are unrelated languages in our sample, we should be careful with any
interpretations here, but let's at least be glad that the algorithm did cluster the Indo-European
languages all together.

.. _KSL.html: examples/KSL.html
