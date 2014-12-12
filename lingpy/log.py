"""Logging utilities"""
from __future__ import unicode_literals, print_function, absolute_import, division
import os
import logging
from logging.config import fileConfig
from tempfile import NamedTemporaryFile
import warnings

from six import text_type

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
level = WARNING
formatter = generic

[formatter_generic]
format = %(asctime)s [%(levelname)s] %(message)s
"""

_logger = None


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
        cfg = Config('logging', default=LOGGING, config_dir=config_dir)
        remove = False
        if cfg.path.exists() and not force_default_config:
            fname = text_type(cfg.path)
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


def info(msg):
    get_logger().info(msg)


def warn(msg):
    get_logger().warn(msg)


def debug(msg):
    get_logger().debug(msg)


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
    logger.warn("Module '{0}' could not be loaded. Some methods may not work properly.".format(name))
