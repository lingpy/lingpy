from __future__ import unicode_literals

from lingpy.tests.util import WithTempDir


class LogTest(WithTempDir):
    def test_new_config(self):
        from lingpy.log import get_logger

        log = get_logger(config_dir=self.tmp, test=True)
        self.assertTrue(hasattr(log, 'info'))

    def test_default_config(self):
        from lingpy.log import get_logger

        log = get_logger(config_dir=self.tmp, force_default_config=True)
        self.assertTrue(hasattr(log, 'info'))

    def test_convenience(self):
        from lingpy.log import info, warn, debug, error, deprecated, missing_module, file_written

        info('m')
        warn('m')
        debug('m')
        error('m')
        deprecated('o', 'n')
        missing_module('m')
        file_written('f')
