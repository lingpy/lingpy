"""
Tests for the read.csv module.
"""
from unittest import TestCase

from six import text_type

from lingpy.compare.lexstat import LexStat
from lingpy.read.csv import *
from lingpy.tests.util import test_data


class Tests(TestCase):
    def setUp(self):
        pass

    def test_read_asjp(self):
        lex = LexStat(read_asjp(
            test_data('asjp_test_list.csv'), family='CELTIC',
            classification='wls_gen'))
        assert len(lex) == 249

        # def evaluate():
        #    return lambda x, y, z: x[y[1]].startswith(z)

        # TODO: Fix value entry, ask Mattis.
        # lex = LexStat(read_asjp(
        #    test_data('asjp_test_list.csv'), family='GERMANIC',
        #    classification='wls_fam,wls_gen', evaluate=evaluate()))

        # assert len(lex) == 1429

        # check if loans have been traced and if at least one word is
        # represented as expected
        # entry = lex.get_dict(doculect='YIDDISH_EASTERN')
        # idx = entry['person'][0]
        # assert lex[idx, 'known_borrowings'] == 1
        # assert lex[idx, 'counterpart'] == "pErzon"

    def test_csv2list(self):
        if_path1 = test_data('test_csv.csv')
        if_path2 = test_data('test_csv')

        # check default setting
        dat1 = csv2list(if_path1)

        # pass data type and header
        dat2 = csv2list(
            if_path2,
            fileformat='csv',
            dtype=[text_type, text_type, text_type, int, text_type],
            sep='\t',
            header=True
        )

        # pass another separator
        dat3 = csv2list(
            if_path1,
            sep='_'
        )

        # modify the comment char
        dat4 = csv2list(
            if_path1,
            comment="?"
        )

        # check for correct parsing
        assert dat1[0][1] == 'is'
        assert dat2[0][1] == 'are'
        assert sum([x[3] for x in dat2]) == 8
        assert dat3[0][0] == 'This\tis\tthe'
        assert dat4[3][0] == '#I'

    def test_csv2dict(self):
        if_path1 = test_data('test_csv.csv')
        if_path2 = test_data('test_csv')

        # check default setting
        dat1 = csv2dict(if_path1)

        # pass data type and header
        dat2 = csv2dict(
            if_path2,
            fileformat='csv',
            dtype=[text_type, text_type, text_type, int, text_type],
            sep='\t',
            header=True
        )

        # pass another separator
        dat3 = csv2dict(
            if_path1,
            sep='_'
        )

        # modify the comment char
        dat4 = csv2dict(
            if_path1,
            comment="?"
        )

        # check for correct results
        assert 'This' in dat1
        assert 'This' not in dat2
        assert dat2['We'][2] == 2
        assert dat3['This\tis\tthe'][0] == 'head\tline'
        assert len(dat4) == 4 and '#I' in dat4

    def test_csv2multidict(self):
        if_path1 = test_data('test_csv.csv')

        md = csv2multidict(if_path1)

        assert md['We']['is'] == 'are'
        assert sum([int(md[x]['head']) for x in md]) == 8
