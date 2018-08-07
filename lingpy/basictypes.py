from functools import partial

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
        super(_strings, self).append(other[0])

    def extend(self, other):
        other = _strings(self._type, other)
        super(_strings, self).extend(other)


strings = partial(_strings, str)
ints = partial(_strings, int)
floats = partial(_strings, float)


class lists(_strings):

    def __init__(self, iterable, sep=" + "):
        _strings.__init__(self, str, iterable)
        self.n = [strings(x) for x in iterable.split(sep)]
        self.sep = sep

    def __add__(self, other):
        return lists(str(self)+self.sep+str(other))

    def extend(self, other):
        super(lists, self).extend(lists('')+other)



        
        
