from __future__ import unicode_literals
from lingpy.basictypes import strings, ints, floats, lists
from nose.tools import assert_raises
from six import text_type

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

    # list check
    assert str(lists('b a + b a') + 'm a') == 'b a + b a + m a'

    # append
    app = strings('1 2 3')
    app.append('4')
    assert str(app) == '1 2 3 4'

    app = ints('1 2 3')
    app.extend('4 5')
    assert str(app) == '1 2 3 4 5'

    assert_raises(ValueError, lambda x: strings(x).append('2 3'), '1 2')

    app = lists('1 2 3 + 1 2')
    app.extend('1 2')
    assert str(app) == '1 2 3 + 1 2 + 1 2'

    app = strings('1 2 3')
    app[1] = 2
    assert app[1] == '2'
