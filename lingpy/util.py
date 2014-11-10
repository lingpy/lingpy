import io
import unicodedata

from pathlib import Path
from six import text_type


def _str_path(path):
    return text_type(path) if isinstance(path, Path) else path


def write_text_file(path, content):
    with io.open(_str_path(path), 'w', encoding='utf8') as fp:
        fp.write(content)


def read_text_file(path, normalize=None, lines=False):
    """Reads a text file encoded in utf-8.

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
