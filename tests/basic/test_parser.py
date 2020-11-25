from collections import defaultdict
from operator import itemgetter
from unittest import TestCase

import pytest

from lingpy.basic.parser import QLCParser, QLCParserWithRowsAndCols
from lingpy.basic.wordlist import Wordlist
from lingpy.cache import path
from lingpy.util import data_path


@pytest.fixture
def parser(test_data):
    return QLCParser(str(test_data / 'KSL.qlc'))


def test_init(test_data):
    p = QLCParser({0: ['a']})
    QLCParser(p)
    with pytest.raises(IOError):
        QLCParser('not-extisting-file')
    with pytest.raises(TypeError):
        QLCParser(None)
    with pytest.raises(ValueError):
        QLCParserWithRowsAndCols({0: ['a']}, 'x', 'y', {})

    with pytest.raises(ValueError):
        QLCParserWithRowsAndCols(
            {0: ['concept', 'language', 'bla'], 1: ['bla', 'blu']}, 'concept', 'language', '')

    p2 = QLCParserWithRowsAndCols(
        str(test_data / 'bad_file2.tsv'), 'concept', 'language',
        data_path('conf', 'wordlist.rc')
    )

    assert p2.get_entries('cogid')[0][-1] == 'ff'
    with pytest.raises(KeyError):
        p2.__getitem__(tuple([2000, 'bla']))
    assert p2[3, 'language'] == 'l3'
    assert p2[3, 'nothing'] is None


def test_getitem(parser):
    key = list(parser._data.keys())[0]
    assert parser[key] == parser._data[key]
    assert parser[key, 'cogid']

    with pytest.raises(KeyError):
        _ = parser[None]
    parser._meta['aaa'] = 3
    assert parser['aaa']

    with pytest.raises(KeyError):
        parser.__getitem__(2000)
    with pytest.raises(KeyError):
        parser.__getitem__(tuple([2000, 'bla']))


def test_get_entries(test_data):
    parser = QLCParserWithRowsAndCols(str(test_data / 'KSL.qlc'), 'concept', 'cogid', {})
    assert parser.get_entries('cogid')


def test_getattr(test_data):
    parser = QLCParserWithRowsAndCols(str(test_data / 'KSL.qlc'), 'concept', 'cogid', {})
    assert parser.cogid


def test_len(parser):
    assert len(parser)


def test_add_entries(parser, mocker):
    with pytest.raises(ValueError):
        parser.add_entries('', 'taxon', lambda t: t.lower())

    parser.add_entries('lTaxon', 'doculect', lambda t: t.lower(), override=True)
    assert 'ltaxon' in parser.entries

    mocker.patch('lingpy.basic.parser.confirm', mocker.Mock(return_value=True))
    parser.add_entries('ltaxon', 'doculect', lambda t: t.lower())

    parser.add_entries('tg', 'doculect,concept',
                            lambda v, id_: ''.join(v))
    with pytest.raises(ValueError):
        parser.add_entries('ll', 'doculect,concept', lambda x, y: 1 / 0)
    parser.add_entries('ti', defaultdict(int), lambda i: i + 1)
    assert 'ti' in parser.entries

    parser.add_entries('tg', defaultdict(int), lambda i: i + 1, override=True)
    parser.add_entries('tg', 'doculect,concept', lambda v, id_: 'abc', override=True)
