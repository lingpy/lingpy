from lingpy.basictypes import strings, ints, floats, lists

class Tests():
    
    string1 = '1 2 3 + 1 2 3'
    string2 = '1 2 3 1 2 3'
    l = lists(string1)
    assert len(l) == 7
    assert len(l.n) == 2
    assert str(l) == string1

    i = ints(string2)
    assert i[0] == 1
    assert str(i) == string2

    f = floats(string2)
    assert float(i[0]) == f[0]
    assert ' '.join([str(fl).split('.')[0] for fl in f]) == string2
    assert str(f).split()[0].startswith('1.')

    s = strings(string1)
    assert str(s) == string1

    # check for types
    s = strings('1 2 3')
    assert str(s + s) == '1 2 3 1 2 3'
    i = ints('1 2 3')
    assert str(i + [1, 2, 3]) == '1 2 3 1 2 3'
