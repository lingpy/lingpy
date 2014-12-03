from __future__ import division
import sys
import io
import unicodedata
from math import ceil
import logging

from pathlib import Path
from six import text_type

from lingpy.log import get_logger


def _str_path(path):
    return text_type(path) if isinstance(path, Path) else path


def write_text_file(path, content, normalize=None):
    """Write a text file encoded in utf-8.

    :param path: File-system path of the file.
    :content: The text content to be written.
    :param normalize: If not `None` a valid unicode normalization mode must be passed.
    """
    log = get_logger()
    with io.open(_str_path(path), 'w', encoding='utf8') as fp:
        fp.write(unicodedata.normalize(normalize, content) if normalize else content)
    log.info("Data has been written to file <{0}>.".format(_str_path(path)))


def read_text_file(path, normalize=None, lines=False):
    """Read a text file encoded in utf-8.

    :param path: File-system path of the file.
    :param normalize: If not `None` a valid unicode normalization mode must be passed.
    :param lines: Flag signalling whether to return a list of lines.
    :return: File content as unicode object or list of lines as unicode objects.

    .. note:: The whole file is read into memory.
    """
    def _normalize(chunk):
        return unicodedata.normalize(normalize, chunk) if normalize else chunk

    with io.open(_str_path(path), 'r', encoding='utf8') as fp:
        if lines:
            return [_normalize(line) for line in fp]
        else:
            return _normalize(fp.read())


def read_config_file(path, **kw):
    """Read lines of a file ignoring commented lines and empty lines. """
    kw['lines'] = True
    lines = [line.strip() for line in read_text_file(path, **kw)]
    return [line for line in lines if line and not line.startswith('#')]


class ProgressBar(object):
    """A progress bar using console output.

    Usage:

    >>> with ProgressBar('here is the title', 50) as pb:
    >>>     for i in range(50):
    >>>         # do stuff
    >>>         pb.update()
    """
    def __init__(self, title, task_count, cols=100):
        self.log = get_logger().getEffectiveLevel() <= logging.INFO
        self.title = title
        self.cols = cols
        self.step = self.cols / task_count
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
