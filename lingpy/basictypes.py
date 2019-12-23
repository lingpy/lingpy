from functools import partial
from six import text_type as str

class _strings(list):

    def __init__(self, type_, iterable):
        
        list.__init__(self, [type_(x) for x in (iterable if not isinstance(iterable, str) else
                iterable.split())])
        self._type = type_

    def __str__(self):
        return ' '.join([str(x) for x in self])

    def __add__(self, other):
        return _strings(self._type, str(self)+' '+str(_strings(self._type, other)))

    def append(self, other):
        other = _strings(self._type, other)
        if len(other) != 1:
            raise ValueError('Use extend() to add more than one item')
        list.append(self, other[0])

    def extend(self, other):
        other = _strings(self._type, other)
        list.extend(self, other)

    def __setitem__(self, index, item):
        list.__setitem__(self, index, self._type(item))

integer = lambda x: int(x) if x else 0
strings = partial(_strings, str)
ints = partial(_strings, int)
floats = partial(_strings, float)

class aligned(_strings):

    def __init__(self, iterable):
        _strings.__init__(self, str, iterable)
    
    @property
    def a(self):
        out = []
        o = False
        for i in self:
            if i == '(':
                o = True
            if not o:
                out += [i]
            if i == ')':
                o = False
        return out

class lists(_strings):

    def __init__(self, iterable, sep=" + "):
        _strings.__init__(self, str, iterable)
        self.n = [strings(x) for x in (' '.join(iterable).split(sep) if not
            isinstance(iterable, str) else iterable.split(sep))]
        self.sep = sep

    def __add__(self, other):
        return lists(str(self)+self.sep+str(other))

    def extend(self, other):
        super(lists, self).extend(lists('')+_strings(str, other))

    def change(self, i, item):
        self.n[i] = _strings(str, item)
        new_s = self.sep.join([' '.join(x) for x in self.n])
        self.__init__(new_s, sep=self.sep)

