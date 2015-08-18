import os
from unittest import TestCase
from operator import itemgetter
from collections import defaultdict

from mock import patch, Mock

from lingpy.tests.util import test_data
from lingpy.basic.parser import QLCParser, QLCParserWithRowsAndCols


class TestParser(TestCase):
    def setUp(self):
        self.parser = QLCParser(test_data('KSL.qlc'))

    def test_init(self):
        p = QLCParser({0: ['a']})
        QLCParser(p)
        self.assertRaises(IOError, QLCParser, 'not-extisting-file')
        self.assertRaises(TypeError, QLCParser, None)
        self.assertRaises(ValueError, QLCParserWithRowsAndCols, {0: ['a']}, 'x', 'y', {})

    def test__tokenize(self):
        self.parser._tokenize(target='ttokens', source='gloss')

    def test_getitem(self):
        key = list(self.parser._data.keys())[0]
        self.assertEquals(self.parser[key], self.parser._data[key])
        assert self.parser[key, 'cogid']
        self.assertRaises(KeyError, itemgetter(None), self.parser)
        self.parser._meta['aaa'] = 3
        assert self.parser['aaa']

    def test_get_entries(self):
        parser = QLCParserWithRowsAndCols(test_data('KSL.qlc'), 'gloss', 'cogid', {})
        assert parser.get_entries('cogid')

    def test_getattr(self):
        parser = QLCParserWithRowsAndCols(test_data('KSL.qlc'), 'gloss', 'cogid', {})
        assert parser.cogid

    def test_cache(self):
        from lingpy.basic.parser import QLCParser
        from lingpy.basic.wordlist import Wordlist
        from lingpy.cache import path

        filename = 'lingpy_test.qlc'
        self.parser.pickle(filename=filename)
        from_cache = QLCParser.unpickle(filename)
        self.assertEqual(self.parser.header, from_cache.header)
        os.remove(str(path(filename)))

        wl = Wordlist(test_data('KSL.qlc'))
        wl.pickle(filename=filename)
        from_cache = Wordlist.unpickle(filename)
        self.assert_(from_cache._class)
        os.remove(str(path(filename)))

    def test_len(self):
        assert len(self.parser)

    def test_add_entries(self):
        self.assertRaises(
            ValueError, self.parser.add_entries, '', 'taxon', lambda t: t.lower())
            
        self.parser.add_entries('lTaxon', 'taxon', lambda t: t.lower(), override=True)
        assert 'ltaxon' in self.parser.entries

        with patch('lingpy.basic.parser.input', Mock(return_value='y')):
            self.parser.add_entries('ltaxon', 'taxon', lambda t: t.lower())

        self.parser.add_entries('tg', 'taxon,gloss', lambda v, id_: ''.join(v))
        self.assertRaises(
            ValueError, self.parser.add_entries, 'll', 'taxon,gloss', lambda x, y: 1 / 0)
        self.parser.add_entries('ti', defaultdict(int), lambda i: i + 1)
        assert 'ti' in self.parser.entries
        self.parser.add_entries('tg', defaultdict(int), lambda i: i + 1, override=True)
        self.parser.add_entries('tg', 'taxon,gloss', lambda v, id_: 'abc', override=True)

    def test__cache(self):
        parser = QLCParserWithRowsAndCols(test_data('KSL.qlc'), 'gloss', 'cogid', {})
        idx = list(parser._data.keys())[0]
        parser._get_cached(idx)
        parser._get_cached(idx)
        parser._clean_cache()
        parser._data.pop(idx)
        self.assertRaises(KeyError, parser._get_cached, idx)
