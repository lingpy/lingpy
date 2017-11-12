===============================
Whats New in Version |version|?
===============================

Version |version| of LingPy does not introduce many new algorithms. Instead it provides the most
stable version of LingPy we had so far. While we are still working on new algorithms and new
approaches we want to introduce at some point, we decided that we want to offer a version that
presents LingPy in a state in which it can be used without creating too much disappointment with
users due to some unforeseen bugs or insufficient testing. For this reason, we have worked hard to
increase test coverage and make sure that the different classes and functions in LingPy run as
smoothly as possible. Apart from potential bug fixes in this 2.6 version, we plan adding new
algorithms to a future 2.7 version to be released some time in 2018. With this version, we hope to
find time and enough assistance of new developers joining the team, so that we can work on a new
version of LingPy in which the complete data structure is refactored, also allowing us to reduce
dependencies like NumPy which bother our users at times.

Bugfixes in Wordlist Class
--------------------------

We fixed some problems related to the ~lingpy.basic.wordlist.Wordlist class in LingPy. First, if you
make the transition from a Wordlist object to a ~lingpy.compare.lexstat.LexStat object, the data
will be deep-copied. As a result, writing::

  >>> wl = Wordlist('myfile.tsv')
  >>> lex = LexStat(wl)
  >>> lex.add_entries('segments', 'tokens', lambda x: x)

will no longer add the same number of new entries to the Wordlist object (as it happened before).
You can now also easily access the header of a wordlist in its current order::

  >>> wl = Wordlist('myfile.tsv')
  >>> wl.columns
  [...]

This is specifically convenient if you construct a Wordlist by creating a dictionary inside a Python
script.

Adding Nexus Output
-------------------

We had nexus output before, but now, Simon Greenhill has helped us to provide a stable export to
both MrBayes and Beast, also including the situation where you want to calculate rates for each
concept class. 

[...]

Orthography Profile Creation
----------------------------

We introduced the commandline of LingPy in version 2.5, but have not really worked on it since then,
as we think that cognate detection analyses should not go on silent by just typing a command into
the terminal. Instead, we encourage users to learn what the major algorithms are, and how they
should be used. We also recommend you to have a look at an `online tutorial <https://github.com/shh-dlce/qmss-2017/tree/master/LingPy>`_ which Johann-Mattis List prepared early in 2017. However, in one particular cases, the commandline actually proved to be very useful, and this is for the creation of orthography profiles (:evobib:`Moran2017`). 

[...]



What's Next?
------------

.. toctree::
   :maxdepth: 1
   
   intro
   examples
   tutorial/index
   docu/index
   reference
   download
