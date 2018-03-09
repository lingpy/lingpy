LingPy
======

|Build Status| |codecov.io| |DOI| |PyPI version| |Documentation|

Authors: `Johann-Mattis List <https://github.com/linguist>`__ (`Max
Planck Institute for the Science of Human
History <http://shh.mpg.de/>`__), `Simon
Greenhill <https://github.com/simongreenhill>`__ (`Max Planck Institute
for the Science of Human History <http://shh.mpg.de/>`__), and `Robert
Forkel <https://github.com/xrotwang>`__ (`Max Planck Institute for the
Science of Human History <http://shh.mpg.de/>`__)

Collaborators: Steven Moran (`Universität
Zürich <http://www.linguistik.uzh.ch/about/mitglieder/moran.html>`__,
`Peter Bouda <http://www.peterbouda.eu/>`__, Johannes Dellert
(`University of
Tübingen <http://www.sfs.uni-tuebingen.de/~gjaeger/evolaemp/index.html>`__),
Taraka Rama (`Centre for Language Technology <http://clt.gu.se/>`__,
Göteborg), and Tiago Tresoldi
(`UFRGS <http://www.ufrgs.br/english/home>`__)

LingPy is a Python Library for Historical Linguistics. It is being
developed for Python 2.7 and Python 3.x using `a single
codebase <https://docs.python.org/3/howto/pyporting.html>`__.

-  All source code is available at: https://github.com/lingpy/lingpy.
-  Documentation can be found at: http://lingpy.org.
-  For a list of papers in which LingPy was applied, see
   `here <https://github.com/lingpy/lingpy/blob/master/PAPERS.md>`__.

Quick Installation
------------------

For our latest stable version, you can simply use pip or easy\_install
for installation:

.. code:: bash

    $ pip install lingpy

or

.. code:: bash

    $ easy_install lingpy

Depending on which easy\_install or pip version you use, either the
Python2 or the Python3 version of LingPy will be installed.

If you want to install LingPy the current GitHub version on your system,
open a terminal and type in the following:

.. code:: bash

    $ git clone https://github.com/lingpy/lingpy/
    $ cd lingpy
    $ python setup.py install

If the last command above returns you some error regarding user
permissions (usually "Errno 13"), you can install LingPy in your home
Python setup:

::

    $ python setup.py install --user

In order to use the library, start an interactive python session and
import LingPy as follows:

.. code:: python

    >>> from lingpy import *

To install LingPy to hack on it, fork the repository on GitHub, open a
terminal and type:

.. code:: bash

    $ git clone https://github.com/<your-github-user>/lingpy/
    $ cd lingpy
    $ python setup.py develop

This will install LingPy in `"development
mode" <http://pythonhosted.org//setuptools/setuptools.html#development-mode>`__,
i.e. you will be able edit the sources in the cloned repository and
import the altered code just as the regular python package.

.. |Build Status| image:: https://travis-ci.org/lingpy/lingpy.svg?branch=master
   :target: https://travis-ci.org/lingpy/lingpy
.. |codecov.io| image:: http://codecov.io/github/lingpy/lingpy/coverage.svg?branch=master
   :target: http://codecov.io/github/lingpy/lingpy?branch=master
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.597082.svg
   :target: https://doi.org/10.5281/zenodo.597082
.. |PyPI version| image:: https://badge.fury.io/py/lingpy.png
   :target: https://badge.fury.io/py/lingpy
.. |Documentation| image:: https://bade.fury.io/py/lingpy.png
   :target: https://lingpy.github.io
