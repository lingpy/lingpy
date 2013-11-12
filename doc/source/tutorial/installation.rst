.. _Installation Instructions:

Installation Instructions
=========================

Basic Installation on Linux and Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend to use the (hopefully) stable version of LingPy (2.0).  In order to install this
version, simply download it, unpack the directory, than "cd" into it, and type in the prompt:

  .. code-block:: bash
  
     $ python setup.py install

Depending on the python version you use, this will either directly install the Python3 version of
LingPy, or Python3 scripts will be automatically converted into a hopefully stable Python2 version
before they are installed. 
If you use Python3 and Cython is also installed on your system, you can install the C-modules
along with the regular LingPy package. Using these modules will result in better performance. In
order to tell the setup-script to install the C-modules, simply type:

  .. code-block:: bash
      
     $ python setup.py install --with-c

You may need sudo-rights to carry out these command.  If the compilation with C-extensions fails,
you may consider using LingPy without C-extensions (it will still work, but in times a bit slower,
since the alternative is written in pure Python). 

Installation on Windows
^^^^^^^^^^^^^^^^^^^^^^^

Our current version of LingPy for Python3 should basically also run on Windows (our tests were
successfull, but it is always possible that we forget some functionalities provided with LingPy). In
order to install LingPy on a Windows machine, we recommend to use the Cygwin terminal and install
LingPy in the same way in which one would otherwise install it on Linux or Mac machines.

Workarounds
^^^^^^^^^^^

To use the library without compiling it, 
you don't have to run the setup-command, but you can directly put it into your Python-path (under
``site-packages``, somewhere in your system).
Alternatively, you can make a symlink in the ``site-packages``-folder, that you call ``lingpy`` 
and link it to ``lingpy-2.0/lingpy/``. For example:

1. Start the Python interpreter (make sure you are using Python3)::
  
  .. code-block:: bash

     $ python

2. At the prompt, locate the site-packages folder:
  
  .. code-block:: python

     >>> import sys
     >>> print(sys.path)
     ['', '/opt/local/Library/Frameworks/Python.framework/Versions/3.2/lib/python3.2', '/opt/local/Library/Frameworks/Python.framework/Versions/3.2/lib/python3.2/site-packages']

3. Create the symlink (you may need sudo):

  .. code-block:: bash
  
     ln -s /path/to/lingpy/lingpy /path/to/site-packages/symlink

4. Test in interpreter:

  .. code-block:: bash

     $ python
  
  .. code-block:: python
  
     >>> import lingpy

Another, simple way, to use LingPy is to include it in your sys-path just before you call the
library::

   >>> import sys
   >>> sys.path.append("path_to_lingpy")


