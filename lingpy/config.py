# *-* coding: utf-8 *-*
"""Configuration management for lingpy.

Various aspects of lingpy can be configured and customized by the user. This is done with
configuration files in the user's config dir.

.. seealso:: https://pypi.python.org/pypi/appdirs/
"""
from __future__ import unicode_literals, print_function, absolute_import, division
import os
import io

from pathlib import Path
from appdirs import user_config_dir
from six import text_type, PY3
from six.moves.configparser import RawConfigParser


DIR = Path(user_config_dir('lingpy'))


class Config(RawConfigParser):
    def __init__(self, name, default=None, **kw):
        """Initialization.

        :param name: Basename for the config file (suffix .ini will be appended).
        :param default: Default content of the config file.
        """
        self.name = name
        self.default = default
        config_dir = kw.pop('config_dir', None) or text_type(DIR)
        RawConfigParser.__init__(self, kw, allow_no_value=True)
        if self.default:
            if PY3:
                fp = io.StringIO(self.default)
            else:
                fp = io.BytesIO(self.default.encode('utf8'))
            self.readfp(fp)

        cfg_path = os.path.join(config_dir, name + '.ini')
        if os.path.exists(cfg_path):
            assert os.path.isfile(cfg_path)
            self.read(cfg_path)
        else:
            if not os.path.exists(config_dir):
                try:
                    os.mkdir(config_dir)
                except OSError:  # pragma: no cover
                    # this happens when run on travis-ci, by a system user.
                    pass
            if os.path.exists(config_dir):
                with open(cfg_path, 'w') as fp:
                    self.write(fp)
        self.path = Path(cfg_path)
