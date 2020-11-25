from lingpy.compare._structure import cv_templates

from lingpy.basic.wordlist import Wordlist


def test_cv_templates(test_data):
    assert cv_templates(Wordlist(str(test_data / 'KSL5.qlc')), 'French', output='markdown')
    patterns, _, sounds = cv_templates(
        Wordlist(str(test_data / 'KSL5.qlc')), 'French', output=None)
    #assert 'V' in [p[0] for p in patterns]
