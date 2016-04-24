from lingpy.sequence.generate import MCPhon, MCBasic


class Tests():
    
    def setUp(self):
        self.words = [
                'hand',
                'fus',
                'kopf',
                'kind',
                'pferd',
                'hund',
                'maus',
                'katze',
                'tier'
                ]

        self.gen = MCPhon(self.words)

    def test_get_string(self):
        string = self.gen.get_string(new=True)
        string2 = self.gen.get_string(new=True, tokens=True)
        
        assert string.replace(' ','') not in self.words
        assert ''.join([x[1] for x in string2]) not in self.words
    
    def test_evaluate_string(self):
        scores = self.gen.evaluate_string('hatze')
        assert scores[1] < 0
        assert scores[0] > 0
