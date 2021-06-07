import pathlib

import pytest


@pytest.fixture
def test_data():
    return pathlib.Path(__file__).parent / 'test_data'
