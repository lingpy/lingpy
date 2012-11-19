QLC-LingPy
======

Available at: https://github.com/lingpy/lingpy

QLC-LingPy is a Python Library for Historical Linguistics that is under development. It is being developed in Python 3.

To use the library in its pre-setup.py stage, git clone the library and either put the library's "lingpy/lingpy" folder your path (or Python path).

Alternatively, you can make a symlink in your Python 3 site-packages folder called "lingpy" and link it to "lingpy/lingpy". For example:

1. Start the Python interpreter (make sure you are using Python 3)

$ python

2. At the prompt, locate the site-packages folder:

>>> import sys
>>> print(sys.path)
['', '/opt/local/Library/Frameworks/Python.framework/Versions/3.2/lib/python3.2', '/opt/local/Library/Frameworks/Python.framework/Versions/3.2/lib/python3.2/site-packages']

3. Create the symlink (you may need sudo):

ln -s /path/to/lingpy/lingpy /path/to/site-packages/symlink

4. Test in interpreter:

$ python
>>> import lingpy














