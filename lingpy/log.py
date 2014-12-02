"""Logging utilities"""
from __future__ import unicode_literals, print_function, absolute_import, division
import os
import logging
from logging.config import fileConfig
from tempfile import NamedTemporaryFile

from six import text_type

from lingpy.config import Config


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
level = WARN
handlers =
qualname = lingpy

[handler_console]
class = StreamHandler
args = (sys.stderr,)
level = NOTSET
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


def file_written(fname, logger=None):
    logger = logger or get_logger()
    logger.info("Data has been written to file <{0}>.".format(fname))
