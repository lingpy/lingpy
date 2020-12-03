import pytest

from lingpy.sequence.generate import MCPhon


@pytest.fixture
def words():
    return [
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


@pytest.fixture
def gen(words):
    return MCPhon(words)


def test_get_string(gen, words):
    string = gen.get_string(new=True)
    string2 = gen.get_string(new=True, tokens=True)

    assert string.replace(' ', '') not in words
    assert ''.join([x[1] for x in string2]) not in words


def test_evaluate_string(gen):
    scores = gen.evaluate_string('hatze')
    assert scores[1] < 0
    assert scores[0] > 0
