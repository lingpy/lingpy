Converting Data to Strings (:py:mod:`~lingpy.convert.strings`)
==============================================================

The strings module provides some general and some specific functions which allow to convert data
into strings which can then be imported by other software tools. You can import it by typing::

   >>> from lingpy.convert.strings import *

Or by typing::
   
   >>> from lingpy.convert import strings

Most of the functions are used internally, being triggered when writing, for example, data from a
~lingpy.basic.wordlist.Wordlist object to file. They can, however, also be used directly, and
especially the ~lingpy.convert.strings.write_nexus function may prove useful to get a more flexible
nexus-output of wordlist data.


Functions
---------
.. autosummary::
   
   ~lingpy.convert.strings.scorer2str
   ~lingpy.convert.strings.msa2str
   ~lingpy.convert.strings.matrix2dst
   ~lingpy.convert.strings.pap2nex
   ~lingpy.convert.strings.pap2csv
   ~lingpy.convert.strings.multistate2nex
   ~lingpy.convert.strings.write_nexus
