# *-* coding utf-8 *-*
from lingpy import *
from collections import defaultdict

class Tests():

    # get all segments
    models = ['sca', 'dolgo', 'art', 'color', 'asjp']
    _model = dict([(model, Model(model)) for model in models])
    segments = set()
    for model in models:
        for segment in _model[model].converter:
            segments.add(segment)
    failures = defaultdict(list)
    for model in models:
        values = list(_model[model].converter.keys())
        for segment in segments:
            if segment not in values and segment != '-':
                failures[model] += [segment]
    if failures:
        error_message = '\n'
        for model in failures:
            error_message += model+'\n'
            error_message += ' // '.join(['"'+x+'"' for x in failures[model]])+'\n\n'
        raise ValueError(error_message)

