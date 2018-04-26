import pytest

from lingpy.algorithm.cython import _talign, _calign, _malign


class TestData:
    def __init__(self):
        self.seq_a = list('abab')
        self.seq_b = list('bababa')
        self.seq_a2 = list('ab1ab')
        self.seq_b2 = list('bab1aba')
        self.m = len(self.seq_a)
        self.n = len(self.seq_b)
        self.scorer = {
            ("a", "b"): -1,
            ("b", "a"): -1,
            ("a", "a"): 1,
            ("b", "b"): 1,
            ("1", "1"): 1,
            ("a", "1"): -1,
            ("b", "1"): -1,
            ("1", "a"): -1,
            ("1", "b"): -1
        }
        self.gop = -1
        self.gap = -1
        self.scale = 0.5
        self.gop_a = [-1, -1, -1, -1]
        self.gop_b = [-2, -1, -1, -1, -1, -0.5]
        self.gop_a2 = [-1, -1, -1, -1, -1]
        self.gop_b2 = [-2, -1, -1, -1, -1, -1, -0.5]
        self.weight_a = [1.5, 1, 1, 1]
        self.weight_b = [1.5, 1, 1, 1, 1, 0.5]
        self.weight_a2 = [1.5, 1, 1, 1, 1]
        self.weight_b2 = [1.5, 1, 1, 1, 1, 1, 0.5]
        self.pro_a = 'abbc'
        self.pro_b = 'abbbbc'
        self.pro_a2 = 'ab1bc'
        self.pro_b2 = 'abb1bbc'

        self.factor = 0.3

        self.profile_a = [['a', 'b'], ['a', 'b']]
        self.profile_b = [['b', 'a', 'b', 'a', 'b', 'a'],
                          ['b', 'a', 'b', 'a', 'X', 'X']]
        self.profile_a2 = [['a', '1', 'b'], ['a', '1', 'b']]
        self.profile_b2 = [['b', 'a', 'b', '1', 'a', 'b', 'a'],
                           ['b', 'a', 'b', '1', 'a', 'X', 'X']]


@pytest.fixture(name='td')
def fixture_test_data():
    return TestData()


def test__talign(td):
    _, _, sim = (_talign.globalign(td.seq_a, td.seq_b, td.m, td.n, td.gop,
                                   td.scale, td.scorer))
    assert round(sim, 2) == 2.5

    _, _, sim = (_talign.semi_globalign(td.seq_a, td.seq_b, td.m, td.n,
                                        td.gop, td.scale, td.scorer))
    assert round(sim, 2) == 4.0

    _, _, sim = (_talign.localign(td.seq_a, td.seq_b, td.m, td.n, td.gop,
                                  td.scale, td.scorer))
    assert round(sim, 2) == 4.0

    _, _, sim = (_talign.dialign(td.seq_a, td.seq_b, td.m, td.n, td.scale,
                                 td.scorer))
    assert round(sim, 2) == 4.0

    _, _, sim = (_talign.globalign(td.seq_b, td.seq_a, td.n, td.m, td.gop,
                                   td.scale, td.scorer))
    assert round(sim, 2) == 2.5

    _, _, sim = (_talign.semi_globalign(td.seq_b, td.seq_a, td.n, td.m, td.gop,
                                        td.scale, td.scorer))
    assert round(sim, 2) == 4.0

    _, _, sim = (_talign.localign(td.seq_b, td.seq_a, td.n, td.m, td.gop,
                                  td.scale, td.scorer))
    assert round(sim, 2) == 4.0

    _, _, sim = (_talign.dialign(td.seq_b, td.seq_a, td.n, td.m, td.scale,
                                 td.scorer))
    assert round(sim, 2) == 4.0

    for mode in ['global', 'overlap', 'local', 'dialign']:
        _, _, sim = (_talign.align_pair(td.seq_a, td.seq_b, td.gop, td.scale,
                                        td.scorer, mode, distance=1))
        assert sim < 1

        _, _, _, dist = (_talign.align_pair(td.seq_a, td.seq_b, td.gop,
                                            td.scale, td.scorer, mode,
                                            distance=2))
        assert dist < 1

        alignments = (_talign.align_pairwise([td.seq_a, td.seq_b], td.gop,
                                             td.scale, td.scorer, mode))
        assert alignments[0][-1] < 1

        alignments = (_talign.align_pairs([[td.seq_a, td.seq_b]], td.gop,
                                          td.scale, td.scorer, mode, 0))
        assert alignments[0][-1] > 1

        alignments = (_talign.align_pairs([[td.seq_a, td.seq_b]], td.gop,
                                          td.scale, td.scorer, mode, 1))
        assert alignments[0][-1] < 1

    for mode in ['global', 'overlap', 'dialign']:
        profiles = (_talign.align_profile(td.profile_a, td.profile_b, td.gop,
                                          td.scale, td.scorer, mode, 0.5))
        assert profiles[-1] == 0.0

    assert _talign.score_profile(['a', 'a'], ['a', 'a'], td.scorer,
                                 td.gop, 0) == 1

    assert _talign.score_profile(['a', 'X'], ['a', 'X'], td.scorer,
                                 td.gop, 0) != 1

    assert _talign.score_profile(['X', 'a'], ['a', 'X'], td.scorer,
                                 td.gop, 0) != 1

    assert _talign.swap_score_profile(['a', '+'], ['X', 'X'], td.scorer,
                                      0, 0) == 0.0

    assert _talign.swap_score_profile(['a', '+'], ['a', 'X'], td.scorer,
                                      0, 0) != 0.0


def test__malign(td):
    _, _, sim = _malign.nw_align(td.seq_a, td.seq_b, td.scorer, td.gap)
    assert sim == 2

    dist = _malign.edit_dist(td.seq_a, td.seq_b, False)
    assert dist == 2

    dist = _malign.edit_dist(td.seq_a, td.seq_b, True)
    assert 1 > dist > 0

    _, _, sim = _malign.sw_align(td.seq_a, td.seq_b, td.scorer, td.gap)
    assert sim == 4

    alms = _malign.we_align(td.seq_a, td.seq_b, td.scorer, td.gap)
    assert alms[0][-1] == 4

    alms = _malign.structalign('abab', 'cdcd')
    assert len(alms[0]) == 1
    assert alms[1] == 2

    sim = _malign.restricted_edit_dist(list('vava'), list('vivi'),
                                       'cvcv', 'cvcv', False)
    assert sim == 2

    sim = _malign.restricted_edit_dist(list('vava'), list('vivi'),
                                       'cvcv', 'cvcv', True)
    assert round(sim, 2) == 0.5


def test__calign(td):
    _, _, sim = (_calign.globalign(td.seq_a, td.seq_b, td.gop_a, td.gop_b,
                                   td.pro_a, td.pro_b, td.m, td.n, td.scale,
                                   td.factor, td.scorer))
    assert round(sim, 2) == 3.1

    _, _, sim = (_calign.semi_globalign(td.seq_a, td.seq_b, td.gop_a, td.gop_b,
                                        td.pro_a, td.pro_b, td.m, td.n,
                                        td.scale, td.factor, td.scorer))
    assert round(sim, 2) == 4.9

    _, _, sim = (_calign.localign(td.seq_a, td.seq_b, td.gop_a, td.gop_b,
                                  td.pro_a, td.pro_b, td.m, td.n, td.scale,
                                  td.factor, td.scorer))
    assert round(sim, 2) == 4.9

    _, _, sim = (_calign.dialign(td.seq_a, td.seq_b, td.pro_a, td.pro_b, td.m,
                                 td.n, td.scale, td.factor, td.scorer))
    assert round(sim, 2) == 4.9

    _, _, sim = (_calign.globalign(td.seq_b, td.seq_a, td.gop_b, td.gop_a,
                                   td.pro_b, td.pro_a, td.n, td.m, td.scale,
                                   td.factor, td.scorer))
    assert round(sim, 2) == 3.1

    _, _, sim = (_calign.semi_globalign(td.seq_b, td.seq_a, td.gop_b, td.gop_a,
                                        td.pro_b, td.pro_a, td.n, td.m,
                                        td.scale, td.factor, td.scorer))
    assert round(sim, 2) == 4.9

    _, _, sim = (_calign.localign(td.seq_b, td.seq_a, td.gop_b, td.gop_a,
                                  td.pro_b, td.pro_a, td.n, td.m, td.scale,
                                  td.factor, td.scorer))
    assert round(sim, 2) == 4.9

    _, _, sim = (_calign.dialign(td.seq_b, td.seq_a, td.pro_b, td.pro_a, td.n,
                                 td.m, td.scale, td.factor, td.scorer))
    assert round(sim, 2) == 4.9

    _, _, sim = (_calign.secondary_globalign(td.seq_a2, td.seq_b2, td.gop_a2,
                                             td.gop_b2, td.pro_a2, td.pro_b2,
                                             len(td.seq_a2), len(td.seq_b2),
                                             td.scale, td.factor, td.scorer,
                                             '1'))
    assert round(sim, 2) == 4.4

    _, _, sim = (_calign.secondary_semi_globalign(td.seq_a2, td.seq_b2,
                                                  td.gop_a2, td.gop_b2,
                                                  td.pro_a2, td.pro_b2,
                                                  len(td.seq_a2),
                                                  len(td.seq_b2), td.scale,
                                                  td.factor, td.scorer, '1'))
    assert round(sim, 2) == 5.2

    _, _, sim = (_calign.secondary_localign(td.seq_a2, td.seq_b2, td.gop_a2,
                                            td.gop_b2, td.pro_a2, td.pro_b2,
                                            len(td.seq_a2), len(td.seq_b2),
                                            td.scale, td.factor,
                                            td.scorer, '1'))
    assert round(sim, 2) == 6.2

    _, _, sim = (_calign.secondary_dialign(td.seq_a2, td.seq_b2, td.pro_a2,
                                           td.pro_b2, len(td.seq_a2),
                                           len(td.seq_b2), td.scale,
                                           td.factor, td.scorer, '1'))
    assert round(sim, 2) == 6.2

    _, _, sim = (_calign.secondary_globalign(td.seq_b2, td.seq_a2, td.gop_b2,
                                             td.gop_a2, td.pro_b2, td.pro_a2,
                                             len(td.seq_b2), len(td.seq_a2),
                                             td.scale, td.factor,
                                             td.scorer, '1'))
    assert round(sim, 2) == 4.4

    _, _, sim = (_calign.secondary_semi_globalign(td.seq_b2, td.seq_a2,
                                                  td.gop_b2, td.gop_a2,
                                                  td.pro_b2, td.pro_a2,
                                                  len(td.seq_b2),
                                                  len(td.seq_a2), td.scale,
                                                  td.factor, td.scorer, '1'))
    assert round(sim, 2) == 5.2

    _, _, sim = (_calign.secondary_localign(td.seq_b2, td.seq_a2, td.gop_b2,
                                            td.gop_a2, td.pro_b2, td.pro_a2,
                                            len(td.seq_b2), len(td.seq_a2),
                                            td.scale, td.factor,
                                            td.scorer, '1'))
    assert round(sim, 2) == 6.2

    _, _, sim = (_calign.secondary_dialign(td.seq_b2, td.seq_a2, td.pro_b2,
                                           td.pro_a2, len(td.seq_b2),
                                           len(td.seq_a2), td.scale, td.factor,
                                           td.scorer, '1'))
    assert round(sim, 2) == 6.2

    for mode in ['global', 'overlap', 'local', 'dialign']:
        _, _, sim = (_calign.align_pair(td.seq_a, td.seq_b, td.weight_a,
                                        td.weight_b, td.pro_a, td.pro_b,
                                        td.gop, td.scale, td.factor, td.scorer,
                                        mode, '1', 1))
        assert sim < 1

        alignments = (_calign.align_pairwise([td.seq_a, td.seq_b],
                                             [td.weight_a, td.weight_b],
                                             [td.pro_a, td.pro_b], td.gop,
                                             td.scale, td.factor, td.scorer,
                                             '1', mode))
        assert alignments[0][-1] < 1

        alignments = (_calign.align_pairwise([td.seq_a2, td.seq_b2],
                                             [td.weight_a2, td.weight_b2],
                                             [td.pro_a2, td.pro_b2], td.gop,
                                             td.scale, td.factor, td.scorer,
                                             '1', mode))
        assert alignments[0][-1] < 1

        alignments = (_calign.align_pairs([[td.seq_a, td.seq_b]],
                                          [[td.weight_a, td.weight_b]],
                                          [[td.pro_a, td.pro_b]], td.gop,
                                          td.scale, td.factor, td.scorer,
                                          mode, '1', 1))
        assert alignments[0][-1] < 1

        alignments = (_calign.align_pairs([[td.seq_a2, td.seq_b2]],
                                          [[td.weight_a2, td.weight_b2]],
                                          [[td.pro_a2, td.pro_b2]], td.gop,
                                          td.scale, td.factor, td.scorer,
                                          mode, '1', 0))
        assert alignments[0][-1] > 1

    for mode in ['global', 'overlap', 'dialign']:
        profiles = (_calign.align_profile(td.profile_a, td.profile_b,
                                          td.weight_a, td.weight_b, td.pro_a,
                                          td.pro_b, td.gop, td.scale, td.factor,
                                          td.scorer, '1', mode, 0.5))
        assert profiles[-1] == 0.0

        profiles = (_calign.align_profile(td.profile_a2, td.profile_b2,
                                          td.weight_a2, td.weight_b2, td.pro_a2,
                                          td.pro_b2, td.gop, td.scale,
                                          td.factor, td.scorer, '1', mode, 0.5))
        assert len(profiles) == 3

    assert _calign.score_profile(['a', 'a'], ['a', 'a'], td.scorer, 0) == 1
    assert _calign.score_profile(['X', 'a'], ['X', 'a'], td.scorer, 0) == 1
    assert _calign.score_profile(['X', 'a'], ['a', 'X'], td.scorer, 0) == 1
    assert _calign.swap_score_profile(['a', '+'], ['X', 'X'], td.scorer, 0) < 0
    assert _calign.swap_score_profile(['a', '+'], ['a', 'X'], td.scorer, 0) < 0


def test_corrdist(td):
    seqs1 = [[td.seq_a, td.seq_b]]
    seqs2 = [[td.seq_a2, td.seq_b2]]
    gops1 = [[td.weight_a, td.weight_b]]
    gops2 = [[td.weight_a2, td.weight_b2]]
    pros1 = [[td.pro_a, td.pro_b]]
    pros2 = [[td.pro_a2, td.pro_b2]]

    for mode in ['global', 'local', 'dialign', 'overlap']:
        corr1 = _calign.corrdist(0.5, seqs1, gops1, pros1, td.gop, td.scale,
                                 td.factor, td.scorer, mode, '1')
        corr2 = _calign.corrdist(0.5, seqs2, gops2, pros2, td.gop,
                                 td.scale, td.factor, td.scorer, mode, '1')
        assert corr1[0]['b', 'b'] == 2
        assert corr2[0]['a', 'a'] == 2
