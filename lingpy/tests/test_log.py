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
