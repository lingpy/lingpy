.. _Installation Instructions:

Installation Instructions
=========================

Since we are still working on the first official release of LingPy-2.0, only a (hopefully) stable development
version is available for download. In order to install this version, simply download it, unpack the directory, 
than "cd" into it, and type in the prompt:

  .. code-block:: bash
  
     $ python setup.py install

You will then be asked whether you want to install LingPy along with the C-modules for better
performance, or not:

  .. code-block:: bash
    
     [i] Do you want to install with C-modules (requires Cython)? (Y/N) 

Just type in 'Y' or 'N', depending on your preferences and the capabilities of your system.
Make sure, you have Python3 installed on your system (the command for calling Python3 may vary). You
may also need sudo-rights to carry out this command. In order to compile the C-modules, you may need Cython. 
If the compilation fails, you may consider using LingPy without C-extensions (it will still work,
but in times a bit slower, since the alternative is written in pure Python). To use the library without compiling it, 
you don't have to run the setup-command, but you can directly put it into your Python-path (under
``site-packages``, somewhere in your system).
Alternatively, you can make a symlink in the ``site-packages``-folder, that you call ``lingpy`` 
and link it to ``lingpy-2.0.dev/lingpy/``. For example:

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


