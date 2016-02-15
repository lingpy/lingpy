from unittest import TestCase
import lingpy.meaning as lpm

class TestBasVoc(TestCase):
    
    def setUp(self):

        self.basvoc = lpm.basvoc.BasVoc()

    def test_init(self):
        basvoc = lpm.basvoc.BasVoc

        self.assertRaises(OSError, basvoc, 'x')
    
    def test_get_list(self):
        jach = self.basvoc.get_list('jachontov','number','item')
        for a,b in [['94', 'water'], ['25', 'eye'], ['45', 'know'], ['86', 'this'], ['84', 'tail'], ['87', 'thou'], ['28', 'fire'], ['89', 'tooth'], ['63', 'one'], ['32', 'full'], ['59', 'new'], ['42', 'I'], ['96', 'what'], ['82', 'sun'], ['61', 'nose'], ['37', 'hand'], ['18', 'dog'], ['24', 'egg'], ['81', 'stone'], ['88', 'tongue'], ['54', 'moon'], ['108', 'wind'], ['98', 'who'], ['104', 'salt'], ['50', 'louse'], ['91', 'two'], ['29', 'fish'], ['21', 'ear'], ['41', 'horn'], ['9', 'blood'], ['17', 'die'], ['110', 'year'], ['57', 'name'], ['10', 'bone'], ['33', 'give']]:

            assert a,b in jach
    def test_get_sublist(self):
        jachdolgo = self.basvoc.get_sublist('jachontov','dolgopolsky','item')
        assert ''.join(sorted(jachdolgo)) == ''.join(sorted(['tooth', 'eye',
            'second person marker', 'tongue', 'two', 'water', 'who/what',
            'who/what', 'first person marker', 'louse', 'name']))


class TestConcepticon(TestCase):
    
    def setUp(self):

        self.cnc = lpm.basvoc.Concepticon()

    def test_compare(self):

        jachdol = self.cnc.compare('jachontov-35','dolgopolsky1964')

        comps = {'1510396': ['louse', 'louse'], '1510400': ['louse', 'louse'], '4118': ['water', 'water'], '5440': ['tongue', 'tongue'], '5450': ['I', 'first person marker'], '5490': ['who', 'who/what'], '5491': ['what', 'who/what'], '5511': ['eye', 'eye'], '5946': ['tooth', 'tooth'], '5995': ['name', 'name'], '6248': ['two', 'two']}

        assert ''.join(sorted([a+b for a,b in jachdol.values()])) == ''.join(sorted([a+b for a,b in
            comps.values()]))

    def test_repr(self):

        assert repr(self.cnc).startswith('<concepticon')

    def test_str(self):

        assert 'dolgopolsky1964' in str(self.cnc)

    def test_getitem(self):

        assert not self.cnc['x']
        assert isinstance(self.cnc['dolgopolsky1964'], dict)
        


