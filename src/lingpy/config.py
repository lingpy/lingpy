"""Configuration management for lingpy.

Various aspects of lingpy can be configured and customized by the user. This is done with
configuration files in the user's config dir.

.. seealso:: https://pypi.python.org/pypi/appdirs/
"""
import io
import configparser
from pathlib import Path

from appdirs import user_config_dir

DIR = Path(user_config_dir('lingpy'))


class Config(configparser.RawConfigParser):
    def __init__(self, name, default=None, **kw):
        """Initialization.

        :param name: Basename for the config file (suffix .ini will be appended).
        :param default: Default content of the config file.
        """
        self.name = name
        self.default = default
        config_dir = Path(kw.pop('config_dir', None) or DIR)
        configparser.RawConfigParser.__init__(self, kw, allow_no_value=True)
        if self.default:
            fp = io.StringIO(self.default)
            self.read_file(fp)

        cfg_path = config_dir.joinpath(name + '.ini')
        if cfg_path.exists():
            assert cfg_path.is_file()
            self.read(str(cfg_path))
        else:
            if not config_dir.exists():
                try:
                    config_dir.mkdir()
                except OSError:  # pragma: no cover
                    # this happens when run on travis-ci, by a system user.
                    pass
            if config_dir.exists():
                with open(cfg_path.as_posix(), 'w') as fp:
                    self.write(fp)
        self.path = cfg_path
