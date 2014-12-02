from __future__ import unicode_literals
import os
from io import open

from lingpy.tests.util import WithTempDir


class ConfigTest(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.cfg = os.path.join(self.tmp, '.config')

    def _make_one(self, **kw):
        from lingpy.config import Config

        return Config('test', config_dir=self.cfg, **kw)

    def test_new_config(self):
        self.assertFalse(os.path.exists(self.cfg))
        cfg = self._make_one()
        self.assertTrue(cfg.path.exists())

    def test_existing_config(self):
        os.mkdir(self.cfg)
        with open(os.path.join(self.cfg, 'test.ini'), 'w', encoding='utf8') as fp:
            fp.write("""\
[section]
option = 12
""")
        cfg = self._make_one()
        self.assertEqual(cfg.get('section', 'option'), '12')

    def test_default(self):
        cfg = self._make_one(default="""\
[section2]
option2 = 7
""")
        self.assertEqual(cfg.get('section2', 'option2'), '7')
