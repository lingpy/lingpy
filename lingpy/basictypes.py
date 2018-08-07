from array import array


class strings(list):

    def __init__(self, iterable):
        
        list.__init__(self, iterable if not isinstance(iterable, str) else
                iterable.split())

    def __str__(self):

        return ' '.join(self)


class ints(list):

    def __init__(self, iterable):

        list.__init__(self, [int(x) for x in (iterable if
           not isinstance(iterable, str) else iterable.split())])

    def __str__(self):

        return ' '.join([str(x) for x in self])


class floats(list):

    def __init__(self, iterable):

        list.__init__(self, [float(x) for x in (iterable if
           not isinstance(iterable, str) else iterable.split())])

    def __str__(self):

        return ' '.join([str(x) for x in self])


class lists(strings):

    def __init__(self, iterable, sep=" + "):

        strings.__init__(self, iterable)
        self.n = [strings(x) for x in iterable.split(sep)]
        self.sep = sep

    def __str__(self):

        return ' '.join(self)
        
        
