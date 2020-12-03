import pytest

from lingpy import Wordlist
from lingpy.evaluate.alr import med


def test_med(test_data):
    wl = Wordlist(str(test_data / 'KSL.qlc'))

    assert med(wl, gold='gloss', test='gloss', classes=False) == pytest.approx(0.0)
    assert med(wl, gold='tokens', test='tokens', classes=True) == pytest.approx(0.0)
