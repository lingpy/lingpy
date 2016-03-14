from lingpy.tests.util import test_data
from lingpy.read.csv import *
from lingpy.meaning.glosses import parse_gloss, compare_conceptlists


def test_parse_gloss():
    data = csv2list(test_data('Ringe-2002-421.tsv'))[1:]
    target_list = [x[0] for x in csv2list(test_data('target_list_ringe.tsv'))]
    for line, target in zip(data, target_list):
        datum = line[1]
        glosses = parse_gloss(datum)
        for a, b, c, d, e, f, g, h, i in glosses:
            # print(datum,'=>',','.join([x for x in [a,b,c,d,e,f,g,''.join(h),i]]),'\t',a)
            assert a == target


def test_conceptlists():
    comparison = compare_conceptlists(test_data('listA.tsv'), test_data('listB.tsv'))
