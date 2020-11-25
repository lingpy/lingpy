"""Logging utilities"""
import os
import sys
import logging
from logging.config import fileConfig
from tempfile import NamedTemporaryFile
import warnings

from .config import Config


LOGGING = """
[loggers]
keys = root, lingpy

[handlers]
keys = console

[formatters]
keys = generic

[logger_root]
level = INFO
handlers = console

[logger_lingpy]
# a level of WARN is equivalent to lingpy's defaults of verbose=False, debug=False
level = INFO
handlers =
qualname = lingpy

[handler_console]
class = StreamHandler
args = (sys.stderr,)
level = INFO
formatter = generic

[formatter_generic]
format = %(asctime)s [%(levelname)s] %(message)s
"""

_logger = None


class CustomFilter(logging.Filter):
    def filter(self, record):
        if getattr(record, 'lines', None):
            for line in record.lines:
                record.msg = record.msg + '\n\t' + '%s' % (line,)
        return super(CustomFilter, self).filter(record)


def get_logger(config_dir=None, force_default_config=False, test=False):
    """Get a logger configured according to the lingpy log config file.

    Note: If no logging configuration file exists, it will be created.

    :param config_dir: Directory in which to look for/create the log config file.
    :param force_default_config: Configure the logger using the default config.
    :param test: Force reconfiguration of the logger.
    :return: A logger.
    """
    global _logger
    if _logger is None or force_default_config or test:
        _logger = logging.getLogger('lingpy')
        _logger.addFilter(CustomFilter())
        testing = len(sys.argv) and sys.argv[0].endswith('nosetests')
        if not (force_default_config or test) and testing:  # pragma: no cover
            _logger.setLevel(logging.CRITICAL)
        else:
            cfg = Config('logging', default=LOGGING, config_dir=config_dir)
            remove = False
            if cfg.path.exists() and not force_default_config:
                fname = str(cfg.path)
            else:
                with NamedTemporaryFile(delete=False) as fp:
                    fp.write(LOGGING.encode('utf8'))
                    fname = fp.name
                    remove = True
            fileConfig(fname, disable_existing_loggers=False)
            if remove:
                os.remove(fname)
    return _logger


def get_level():
    return get_logger().getEffectiveLevel()


def info(msg, **kw):
    get_logger().info(msg, **kw)


def warning(msg):
    get_logger().warning(msg)


def debug(msg, **kw):
    get_logger().debug(msg, **kw)


def error(msg, **kw):
    get_logger().error(msg, **kw)


def file_written(fname, logger=None):
    logger = logger or get_logger()
    logger.info("Data has been written to file <{0}>.".format(fname))


def deprecated(old, new):
    warnings.warn(
        "Use of '{0}' is deprecated, use '{1}' instead.".format(old, new),
        DeprecationWarning)


def missing_module(name, logger=None):
    logger = logger or get_logger()
    logger.warning(
        "Module '{0}' could not be loaded. Some methods may not work properly.".format(
            name))


class Logging(object):
    """
    A context manager to execute a block of code at a specific logging level.
    """
    def __init__(self, level=logging.DEBUG, logger=None):
        self.level = level
        self.logger = logger or get_logger()
        self.linpy_logger_level = self.logger.getEffectiveLevel()
        root = logging.getLogger()
        self.root_logger_level = root.getEffectiveLevel()
        self.root_handler_level = root.handlers[0].level

    def __enter__(self):
        self.logger.setLevel(self.level)
        root = logging.getLogger()
        root.setLevel(self.level)
        root.handlers[0].setLevel(self.level)
        return self.logger

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.setLevel(self.linpy_logger_level)
        root = logging.getLogger()
        root.setLevel(self.root_logger_level)
        root.handlers[0].setLevel(self.root_handler_level)
