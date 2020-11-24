import pytest

from lingpy import *


@pytest.mark.parametrize(
    'model',
    ['sca', 'dolgo', 'art', 'color', 'asjp', 'cv']
)
def test_model(model):
    # get all segments
    _model = Model(model)

    values = list(_model.converter.keys())
    for segment in set(seg for seg in _model.converter):
        assert (segment in values) or (segment == '-')
