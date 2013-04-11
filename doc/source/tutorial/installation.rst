.. _Installation Instructions:

Installation Instructions
=========================

Since we are still working on the first official release of LingPy-2.0 for Python3, an automatic installation
with a setup.py or tools like easy_install is currently not supported.  To use the library in this
pre-setup.py stage, git clone the library from https://github.com/lingpy/lingpy and put the
library's ``lingpy/lingpy`` folder in your path (or Python path).  Alternatively, you can make a
symlink in your Python3 ``site-packages`` folder called ``lingpy`` and link it to ``lingpy/lingpy``. For
example:

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


