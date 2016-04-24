"""
LingPy comes along with many different kinds of predefined data. When loading
the library, the following dictionary is automatically loaded and employed by
all LingPy modules:

    .. py:data:: rcParams : dict

       As an alternative to all global variables, this dictionary
       contains all these variables, and additional ones. This dictionary is used
       for internal coding purposes and stores parameters that are globally set (if
       not defined otherwise by the user), such as
          
          * specific debugging messages (warnings, messages, errors)
          * default values, such as "gop" (gap opening penalty), "scale" (scaling factor
          * by which extended gaps are penalized), or "figsize" (the default size of
          * figures if data is plotted using matplotlib).
          
       These default values can be changed with help of the ``rc`` function that takes any
       keyword and any variable as input and adds or modifies the specific key of the
       rcParams dictionary, but also provides more complex functions that change whole
       sets of variables, such as the following statement::
          
          >>> rc(schema="asjp")
          
       which switches the variables "asjp", "dolgo", etc. to the ASCII-based
       transcription system of the ASJP project.

       If you want to change the content of c{rcParams} directly, you need to
       import the dictionary explicitly::

          >>> from lingpy.settings import rcParams

       However, changing the values in the dictionary randomly can produce
       unexpected behavior and we recommend to use the regular ``rc`` function
       for this purpose.


.. autofunction:: lingpy.settings.rc

"""

from .model import *
