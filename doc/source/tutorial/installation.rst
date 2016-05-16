.. _Installation Instructions:

Installation Instructions
=========================

Basic Installation on Linux and Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend to use the stable version of LingPy |version|. In order to install this
version, you can use PIP:

  .. code-block:: bash

     $ pip install lingpy

Alternatively, simply download the source code, unpack the directory, than "cd" into it, and type in the prompt:

  .. code-block:: bash
  
     $ python setup.py install
 
If you use Python3 and Cython is also installed on your system, you can install the C-modules
along with the regular LingPy package. Using these modules will result in better performance. In
order to tell the setup-script to install the C-modules, simply type:

  .. code-block:: bash
      
     $ python setup.py install --with-c

You may need sudo-rights to carry out these command.  If the compilation with C-extensions fails,
you may consider using LingPy without C-extensions (it will still work, but in times a bit slower,
since the alternative is written in pure Python). 


