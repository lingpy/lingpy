import pathlib

import pytest


@pytest.fixture
def tmppath(tmpdir):
    return pathlib.Path(str(tmpdir))


@pytest.fixture
def test_data():
    return pathlib.Path(__file__).parent / 'test_data'
