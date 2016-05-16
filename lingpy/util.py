from __future__ import division, unicode_literals
import sys
import io
import unicodedata
from math import ceil
import logging
import os
from tempfile import NamedTemporaryFile
from functools import partial
import json
import itertools
import types
import csv

from pathlib import Path
from six import text_type, PY3
from six.moves import input

import lingpy
from lingpy.log import get_level, file_written
from lingpy.settings import rcParams


PROG = 'LingPy-{0}'.format(lingpy.__version__)


def charstring(id_, char='X', cls='-'):
    return '{0}.{1}.{2}'.format(id_, char, cls)


def combinations2(iterable):
    """
    Convenience shortcut
    """
    return itertools.combinations(iterable, 2)

def multicombinations2(iterable):
    """
    Convenience shortcut, for the name, see the Wikipedia article on Combination.

    https://en.wikipedia.org/wiki/Combination#Number_of_combinations_with_repetition
    """
    return itertools.combinations_with_replacement(iterable, 2)


def join(sep, *args, **kw):
    """
    Convenience shortcut. Strings to be joined do not have to be passed as list or tuple.
    
    Notes
    -----
    An implicit conversion of objects to strings is performed as well.
    
    """
    if len(args) == 1 and isinstance(args[0], (list, tuple, types.GeneratorType)):
        args = args[0]
    condition = kw.get('condition', lambda x: True)
    return sep.join(['%s' % arg for arg in args if condition(arg)])


dotjoin = partial(join, '.')
tabjoin = partial(join, '\t')

def confirm(question, default=False):
    """
    Ask user a yes/no question and return their response as True or False.
    
    Parameters
    ----------
    question : str
        The question you want to answer.

    Notes
    -----
    ``question`` should be a simple, grammatically complete question such as
    "Do you wish to continue?", and will have a string similar to " [Y/n] "
    appended automatically. This function will *not* append a question mark for
    you.

    By default, when the user presses Enter without typing anything, or the user
    is not recognized as "yes" or "no", "no" is assumed. This can be changed by
    specifying ``default=True``.

    Adapted from fabric.contrib.console.confirm
    """
    response = input("%s [%s] " % (question, "Y/n" if default else "y/N")).lower()
    if response in rcParams['answer_yes']:
        return True

    if response in ['n', 'no']:
        return False

    # Didn't get empty, yes or no, so complain and loop
    return default


class TemporaryPath(object):
    def __init__(self, suffix=''):
        fp = NamedTemporaryFile(suffix=suffix)
        self.name = fp.name
        fp.close()

    def __enter__(self):
        return self.name

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.name):
            os.remove(self.name)


def lingpy_path(*comps):
    return os.path.join(os.path.dirname(lingpy.__file__), *comps)

data_path = partial(lingpy_path, 'data')


def _str_path(path, mkdir=False):
    """Get a file-system path as text_type, suitable for passing into io.open.
    
    Parameters
    ----------
    path : {text_type, Path}
        A fs path either as Path instance or as text_type.
    mkdir : bool (default=False)
        If True, create the directories within the path.
    
    Returns
    -------
    path : text_type
        The path as text_type.
    """
    res = text_type(path) if isinstance(path, Path) else path
    if mkdir:
        dirname = os.path.dirname(res)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)
    return res


def write_text_file(path, content, normalize=None, log=True):
    """Write a text file encoded in utf-8.
    
    Parameters
    ----------
    path : str
        File-system path of the file.
    content : str
        The text content to be written.
    normalize : { None, "NFC", "NFD" } (default=False)
        If not `None` a valid unicode normalization mode must be passed.
    log : bool (default=True)
        Indicate whether you want to log the result of the file writing
        process.

    """
    if not isinstance(content, text_type):
        content = lines_to_text(content)
    with io.open(_str_path(path, mkdir=True), 'w', encoding='utf8') as fp:
        fp.write(unicodedata.normalize(normalize, content) if normalize else content)
    if log:
        file_written(_str_path(path))

def lines_to_text(lines):
    return ''.join(line if line.endswith('\n') else line + '\n' for line in lines)


class TextFile(object):
    def __init__(self, path, log=True):
        self.path = path
        self.log = log
        self.fp = io.open(_str_path(path, mkdir=True), 'w', encoding='utf8')

    def __enter__(self):
        return self.fp

    def __exit__(self, type, value, traceback):
        self.fp.close()
        if self.log:
            file_written(_str_path(self.path))

def read_csv_file(path, delimiter=",", quotechar='"',
        normalize=None):
    """
    Load a normal CSV file.

    Parameters
    ----------
    path : str
        The path to your CSV file.
    delimiter : str
        The delimiter in the CSV file.
    quotechar : str
        The quote character in your data.
    normalize : { None, "NFC", "NFC" }
        If not `None` a valid unicode normalization mode must be passed.

    """
    def _normalize(chunk):
        return chunk #unicodedata.normalize(normalize, chunk) if normalize else chunk
    
    if PY3:
        with io.open(_str_path(path), 'r', encoding='utf8') as fp:
            reader = csv.reader(fp, delimiter=delimiter,
                quotechar=quotechar)
            return [[_normalize(cell) for cell in row] for row in reader]
    else:
        with open(_str_path(path), 'rb') as fp:
            reader = csv.reader(fp, delimiter=bytes(delimiter),
                    quotechar=bytes(quotechar))
            return [row for row in reader]



def read_text_file(path, normalize=None, lines=False):
    """
    Read a text file encoded in utf-8.

    Parameters
    ----------
    path : { Path, str }
        File-system path of the file.
    normalize : { None, "NFC", "NFC" }
        If not `None` a valid unicode normalization mode must be passed.
    lines : bool (default=False)
        Flag signalling whether to return a list of lines (without
        the line-separation character).
    
    Returns
    -------
    file_content : { list, str }
        File content as unicode object or list of lines as unicode objects.
    
    Notes
    -----
    The whole file is read into memory.
    
    """
    def _normalize(chunk):
        return unicodedata.normalize(normalize, chunk) if normalize else chunk

    with io.open(_str_path(path), 'r', encoding='utf8') as fp:
        if lines:
            return [_normalize(line.strip('\r\n')) for line in fp]
        else:
            return _normalize(fp.read())

def as_string(obj, pprint=False):
    obj = text_type(obj)
    if get_level() <= logging.INFO and pprint:
        print(obj)
    return obj

def read_config_file(path, **kw):
    """Read lines of a file ignoring commented lines and empty lines. """
    kw['lines'] = True
    lines = [line.strip() for line in read_text_file(path, **kw)]
    return [line for line in lines if line and not line.startswith('#')]


class ProgressBar(object):
    """A progress bar using console output.
    
    Notes
    -----
    Usage example::

        >>> with ProgressBar('here is the title', 50) as pb:
        >>>     for i in range(50):
        >>>         # do stuff
        >>>         pb.update()

    """
    def __init__(self, title, task_count, cols=100):
        self.log = get_level() <= logging.INFO
        self.title = title
        self.cols = cols
        self.step = cols if not task_count else self.cols / task_count
        # number of columns we have already written:
        self.written = 0
        # number of tasks
        self.count = 0

    def __enter__(self):
        if self.log:
            sys.stdout.write(
                '|' + ' {0} '.format(self.title).center(self.cols, '-') + '|\r|')
        return self

    def __exit__(self, type, value, traceback):
        if self.log:
            sys.stdout.write('|\n')

    def update(self):
        self.count += 1
        # compute how many of the columns should have been written by now:
        percentage = int(ceil(self.count * self.step))
        if percentage > self.written:
            # and write what's missing:
            to_write = '+' * (percentage - self.written)
            if self.log and self.written < self.cols:
                sys.stdout.write(to_write)
                sys.stdout.flush()
            self.written += len(to_write)


def setdefaults(d, **kw):
    """Shortcut for a common idiom, setting multiple default values at once.
    
    Parameters
    ----------
    d : dict
        Dictionary to be updated.
    kw : dict
        Dictionary with default values.
    """
    for k, v in kw.items():
        d.setdefault(k, v)


class cached_property(object):

    """
    Decorator for read-only properties evaluated only once.
    
    Notes
    -----
    It can be used to create a cached property like this::

        import random

        # the class containing the property must be a new-style class
        class MyClass(object):
            # create property whose value is cached
            @cached_property()
            def randint(self):
                # will only be evaluated once.
                return random.randint(0, 100)

    The value is cached  in the '_cache' attribute of the object instance that
    has the property getter method wrapped by this decorator. The '_cache'
    attribute value is a dictionary which has a key for every property of the
    object which is wrapped by this decorator. Each entry in the cache is
    created only when the property is accessed for the first time and is the last
    computed property value.

    To expire a cached property value manually just do::

        del instance._cache[<property name>]

    Inspired by the recipe by Christopher Arndt in the PythonDecoratorLibrary
    """

    def __call__(self, fget):
        self.fget = fget
        self.__doc__ = fget.__doc__
        self.__name__ = fget.__name__
        self.__module__ = fget.__module__
        return self

    def __get__(self, inst, owner):
        if inst is None:
            return
        if not hasattr(inst, '_cache'):
            inst._cache = {}
        if self.__name__ not in inst._cache:
            inst._cache[self.__name__] = self.fget(inst)
        return inst._cache[self.__name__]


def identity(x):
    return x


def jsondump(obj, path, **kw):
    """python 2 + 3 compatible version of json.dump.
    
    Parameters
    ----------
    obj : object
        The object to be dumped.
    path : { str, Path }
        The path of the JSON file to be written.
    """
    _kw = dict(mode='w')
    if PY3:  # pragma: no cover
        _kw['encoding'] = 'utf8'
    with open(path, **_kw) as fp:
        return json.dump(obj, fp, **kw)


def jsonload(path, **kw):
    """python 2 + 3 compatible version of json.load.
    
    Returns
    -------
    obj : object
        The python object read from path.
    """
    _kw = {}
    if PY3:  # pragma: no cover
        _kw['encoding'] = 'utf8'
    with open(path, **_kw) as fp:
        return json.load(fp, **kw)
