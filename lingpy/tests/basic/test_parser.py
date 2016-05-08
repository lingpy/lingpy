import os
from unittest import TestCase
from operator import itemgetter
from collections import defaultdict

from mock import patch, Mock

from lingpy.tests.util import test_data
from lingpy.util import data_path
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
        
        self.assertRaises(ValueError, QLCParserWithRowsAndCols, 
                {0: ['concept', 'language', 'bla'],
                    1 : ['bla', 'blu']}, 'concept', 'language', '')
        
        p2 = QLCParserWithRowsAndCols(test_data('bad_file2.tsv'), 'concept',
            'language', data_path('conf', 'wordlist.rc'))
        assert p2.get_entries('cogid')[0][-1] == 'ff'
        self.assertRaises(KeyError, p2.__getitem__, tuple([2000, 'bla']))
        assert p2[3, 'language'] == 'l3'
        assert p2[3, 'nothing'] is None

    def test_getitem(self):
        key = list(self.parser._data.keys())[0]
        self.assertEquals(self.parser[key], self.parser._data[key])
        assert self.parser[key, 'cogid']
        self.assertRaises(KeyError, itemgetter(None), self.parser)
        self.parser._meta['aaa'] = 3
        assert self.parser['aaa']
        self.assertRaises(KeyError, self.parser.__getitem__, 2000)
        self.assertRaises(KeyError, self.parser.__getitem__, tuple([2000,
            'bla']))

    def test_get_entries(self):
        parser = QLCParserWithRowsAndCols(test_data('KSL.qlc'), 'concept', 'cogid', {})
        assert parser.get_entries('cogid')

    def test_getattr(self):
        parser = QLCParserWithRowsAndCols(test_data('KSL.qlc'), 'concept', 'cogid', {})
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
            
        self.parser.add_entries('lTaxon', 'doculect', lambda t: t.lower(), override=True)
        assert 'ltaxon' in self.parser.entries

        with patch('lingpy.basic.parser.confirm', Mock(return_value=True)):
            self.parser.add_entries('ltaxon', 'doculect', lambda t: t.lower())

        self.parser.add_entries('tg', 'doculect,concept', lambda v, id_: ''.join(v))
        self.assertRaises(
            ValueError, self.parser.add_entries, 'll', 'doculect,concept', lambda x, y: 1 / 0)
        self.parser.add_entries('ti', defaultdict(int), lambda i: i + 1)
        assert 'ti' in self.parser.entries
        self.parser.add_entries('tg', defaultdict(int), lambda i: i + 1, override=True)
        self.parser.add_entries('tg', 'doculect,concept', lambda v, id_: 'abc', override=True)

