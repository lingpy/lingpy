from __future__ import unicode_literals
import os
from unittest import TestCase


class TestCache(TestCase):
    def test_cache(self):
        from lingpy import cache

        d = {'a': 123}
        filename = 'lingpy_test.CSV'
        cache.dump(d, filename)
        self.assertTrue(cache.path(filename).exists())
        self.assertEqual(cache.load(filename), d)
        os.remove(str(cache.path(filename)))
