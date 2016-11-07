from __future__ import division, unicode_literals
import sys
import io
import unicodedata
from math import ceil
import logging
from tempfile import NamedTemporaryFile
from functools import partial
import itertools
import types

from six import text_type
from clldutils.path import Path, remove, path_component
from clldutils import clilib

import lingpy
from lingpy.log import get_level, file_written


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
confirm = partial(clilib.confirm, default=False)


class TemporaryPath(object):
    def __init__(self, suffix=''):
        fp = NamedTemporaryFile(suffix=suffix)
        self.name = Path(fp.name)
        fp.close()

    def __enter__(self):
        return self.name.as_posix()

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.name.exists():
            remove(self.name)


def lingpy_path(*comps):
    return Path(lingpy.__file__).parent.joinpath(*comps).as_posix()

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
    res = Path(path_component(path))
    if mkdir and res.parent and not res.parent.exists():
        res.parent.mkdir(parents=True)
    return res.as_posix()


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


def identity(x):
    return x
