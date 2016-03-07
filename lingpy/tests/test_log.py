from __future__ import unicode_literals

from lingpy.tests.util import WithTempDir

from lingpy import log


class LogTest(WithTempDir):
    def tearDown(self):
        # Reset the global _logger to make sure, it's recreated with appropriate config
        # when the other tests are run.
        log._logger = None
        WithTempDir.tearDown(self)

    def test_new_config(self):
        l = log.get_logger(config_dir=self.tmp, test=True)
        self.assertTrue(hasattr(l, 'info'))

    def test_default_config(self):
        l = log.get_logger(config_dir=self.tmp, force_default_config=True)
        self.assertTrue(hasattr(l, 'info'))

    def test_convenience(self):
        from lingpy.log import info, warn, debug, error, deprecated, missing_module, file_written

        info('m')
        warn('m')
        debug('m')
        error('m')
        deprecated('o', 'n')
        missing_module('m')
        file_written('f')
