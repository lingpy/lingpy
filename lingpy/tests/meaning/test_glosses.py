from lingpy.tests.util import test_data
from lingpy.read.csv import *
from lingpy.meaning.glosses import parse_gloss

def test_parse_gloss():

    data = csv2list(test_data('Ringe-2002-421.tsv'))[1:]
    target_list = [x[0] for x in csv2list(test_data('target_list_ringe.tsv'))]
    for line,target in zip(data,target_list):
        datum = line[1]
        glosses = parse_gloss(datum)
        for a,b,c,d,e,f in glosses:
            print(datum,'=>',','.join([x for x in [a,b,c,d,e,f]]),'\t',a)
            assert a == target

if __name__ == '__main__':
    test_parse_gloss()

