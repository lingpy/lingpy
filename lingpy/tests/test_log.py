from __future__ import unicode_literals

import logging

from lingpy import log
from lingpy.log import info, warning, debug, error, deprecated, missing_module, \
    file_written
from lingpy.tests.util_testing import WithTempDir, capture_all


class LogTest(WithTempDir):
    def tearDown(self):
        # Reset the global _logger to make sure, it's recreated with appropriate
        # config when the other tests are run.
        log._logger = None
        WithTempDir.tearDown(self)

    def test_Logging_context_manager(self):
        # Note: We must make sure to acquire a fresh logger *within* the capture
        # all context, because otherwise, stderr redirection won't work. Thus,
        # we have to initialize the logger and call one of its logging methods
        # within one function which we will then pass to capture_all.
        def _log(method='debug', with_context_manager=False,
                 level=logging.DEBUG):
            logger = log.get_logger(test=True, force_default_config=True,
                                    config_dir=self.tmp_path())

            method = getattr(logger, method)

            if with_context_manager:
                with log.Logging(logger=logger, level=level):
                    method('')
            else:
                method('')

        with capture_all(_log, with_context_manager=False) as res:
            self.assertNotIn('DEBUG', res[2])

        with capture_all(_log, with_context_manager=True) as res:
            self.assertIn('DEBUG', res[2])

        with capture_all(_log, with_context_manager=True,
                         level=logging.WARN) as res:
            self.assertNotIn('DEBUG', res[2])

        with capture_all(
            _log, method='warning', with_context_manager=True, level=logging.WARN
        ) as res:
            self.assertIn('WARN', res[2])

    def test_new_config(self):
        new_cfg = log.get_logger(config_dir=self.tmp.as_posix(), test=True)
        self.assertTrue(hasattr(new_cfg, 'info'))

    def test_default_config(self):
        default_cfg = log.get_logger(config_dir=self.tmp.as_posix(),
                                     force_default_config=True)
        self.assertTrue(hasattr(default_cfg, 'info'))

    @staticmethod
    def test_convenience():
        info('m')
        warning('m')
        debug('m')
        error('m')
        deprecated('o', 'n')
        missing_module('m')
        file_written('f')
