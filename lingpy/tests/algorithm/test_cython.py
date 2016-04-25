from lingpy.algorithm.cython import _talign
from lingpy.algorithm.cython import _calign
from lingpy.algorithm.cython import _malign

class Tests():

    def setUp(self):

        self.seqA = list('abab')
        self.seqB = list('bababa')
        self.seqA2 = list('ab1ab')
        self.seqB2 = list('bab1aba')
        self.m = len(self.seqA)
        self.n = len(self.seqB)
        self.scorer = {
                ("a", "b") : -1,
                ("b", "a") : -1,
                ("a", "a") : 1,
                ("b", "b") : 1,
                ("1", "1") : 1,
                ("a", "1") : -1,
                ("b", "1") : -1,
                ("1", "a") : -1,
                ("1", "b") : -1
                }
        self.gop = -1
        self.gap = -1
        self.scale = 0.5
        self.gopA = [-1, -1, -1, -1]
        self.gopB = [-2, -1, -1, -1, -1, -0.5]
        self.gopA2 = [-1, -1, -1, -1, -1]
        self.gopB2 = [-2, -1, -1, -1, -1, -1, -0.5]
        self.weightA = [1.5, 1, 1, 1]
        self.weightB = [1.5, 1, 1, 1, 1, 0.5]
        self.weightA2 = [1.5, 1, 1, 1, 1]
        self.weightB2 = [1.5, 1, 1, 1, 1, 1, 0.5]
        self.proA = 'abbc'
        self.proB = 'abbbbc'
        self.proA2 = 'ab1bc'
        self.proB2 = 'abb1bbc'

        self.factor = 0.3

        self.profileA = [['a', 'b'], ['a', 'b']]
        self.profileB = [['b', 'a', 'b', 'a', 'b', 'a'], 
                ['b', 'a', 'b', 'a', 'X', 'X']]
        self.profileA2 = [['a', '1', 'b'], ['a', '1', 'b']]
        self.profileB2 = [['b', 'a', 'b', '1', 'a', 'b', 'a'], 
                ['b', 'a', 'b', '1', 'a', 'X', 'X']]

    def test__talign(self):

        almA, almB, sim = _talign.globalign(self.seqA, self.seqB, self.m,
                self.n, self.gop, self.scale, self.scorer)
        assert round(sim, 2) == 2.5
        almA, almB, sim = _talign.semi_globalign(self.seqA, self.seqB, self.m,
                self.n, self.gop, self.scale, self.scorer)
        assert round(sim, 2) == 4.0
        almA, almB, sim = _talign.localign(self.seqA, self.seqB, self.m,
                self.n, self.gop, self.scale, self.scorer)
        assert round(sim, 2) == 4.0
        almA, almB, sim = _talign.dialign(self.seqA, self.seqB, self.m,
                self.n, self.scale, self.scorer)
        assert round(sim, 2) == 4.0
        
        for mode in ['global', 'overlap', 'local', 'dialign']:
            almA, almB, sim = _talign.align_pair(self.seqA, self.seqB, 
                    self.gop, self.scale, self.scorer, mode, distance=1)
            assert sim < 1
            
            alignments = _talign.align_pairwise([self.seqA, self.seqB],
                    self.gop, self.scale, self.scorer, mode)
            assert alignments[0][-1] < 1
            alignments = _talign.align_pairs([[self.seqA, self.seqB]],
                    self.gop, self.scale, self.scorer, mode, 0)
            assert alignments[0][-1] > 1
            alignments = _talign.align_pairs([[self.seqA, self.seqB]],
                    self.gop, self.scale, self.scorer, mode, 1)
            assert alignments[0][-1] < 1
        for mode in ['global', 'overlap', 'dialign']:
            profiles = _talign.align_profile(self.profileA, self.profileB,
                    self.gop, self.scale, self.scorer, mode, 0.5)
            assert profiles[-1] == 0.0

        assert _talign.score_profile(
                ['a', 'a'], ['a', 'a'], self.scorer, self.gop, 0
                ) == 1

        assert _talign.swap_score_profile(
                ['a', '+'], ['X', 'X'], self.scorer, 0, 0
                ) == 0.0

    def test__malign(self):

        almA, almB, sim = _malign.nw_align(self.seqA, self.seqB, self.scorer, self.gap)
        assert sim == 2
        dist = _malign.edit_dist(self.seqA, self.seqB, False)
        assert dist == 2
        dist = _malign.edit_dist(self.seqA, self.seqB, True)
        assert dist < 1 and dist > 0
        almA, almB, sim = _malign.sw_align(self.seqA, self.seqB, self.scorer, self.gap)
        assert sim == 4
        alms = _malign.we_align(self.seqA, self.seqB, self.scorer, self.gap)
        assert alms[0][-1] == 4
        alms = _malign.structalign('abab', 'cdcd')
        assert len(alms[0]) == 1
        assert alms[1] == 2
        sim = _malign.restricted_edit_dist(list('vava'),
                list('vivi'), 'cvcv', 'cvcv', False)
        assert sim == 2
        sim = _malign.restricted_edit_dist(list('vava'),
                list('vivi'), 'cvcv', 'cvcv', True)
        assert round(sim, 2) == 0.5

    def test__calign(self):

        almA, almB, sim = _calign.globalign(self.seqA, self.seqB, self.gopA,
                self.gopB, self.proA, self.proB, self.m,
                self.n, self.scale, self.factor, self.scorer)
        assert round(sim, 2) == 3.1
        almA, almB, sim = _calign.semi_globalign(self.seqA, self.seqB, self.gopA,
                self.gopB, self.proA, self.proB, self.m,
                self.n, self.scale, self.factor, self.scorer)
        assert round(sim, 2) == 4.9
        almA, almB, sim = _calign.localign(self.seqA, self.seqB, self.gopA,
                self.gopB, self.proA, self.proB, self.m,
                self.n, self.scale, self.factor, self.scorer)
        assert round(sim, 2) == 4.9
        almA, almB, sim = _calign.dialign(self.seqA, self.seqB,
                self.proA, self.proB, self.m,
                self.n, self.scale, self.factor, self.scorer)
        assert round(sim, 2) == 4.9

        almA, almB, sim = _calign.secondary_globalign(self.seqA2, self.seqB2, self.gopA2,
                self.gopB2, self.proA2, self.proB2, len(self.seqA2),
                len(self.seqB2), self.scale, self.factor, self.scorer, '1')
        assert round(sim, 2) == 4.4
        almA, almB, sim = _calign.secondary_semi_globalign(self.seqA2, self.seqB2, self.gopA2,
                self.gopB2, self.proA2, self.proB2, len(self.seqA2),
                len(self.seqB2), self.scale, self.factor, self.scorer, '1')
        assert round(sim, 2) == 5.2
        almA, almB, sim = _calign.secondary_localign(self.seqA2, self.seqB2, self.gopA2,
                self.gopB2, self.proA2, self.proB2, len(self.seqA2),
                len(self.seqB2), self.scale, self.factor, self.scorer, '1')
        assert round(sim, 2) == 6.2
        almA, almB, sim = _calign.secondary_dialign(self.seqA2, self.seqB2,
                self.proA2, self.proB2, len(self.seqA2),
                len(self.seqB2), self.scale, self.factor, self.scorer, '1')
        assert round(sim, 2) == 6.2

        for mode in ['global', 'overlap', 'local', 'dialign']:
            almA, almB, sim = _calign.align_pair(self.seqA, self.seqB, 
                    self.weightA, self.weightB, self.proA, self.proB, self.gop,
                    self.scale, self.factor, self.scorer, mode,
                    '1', 1)
            assert sim < 1
            
            alignments = _calign.align_pairwise([self.seqA, self.seqB],
                    [self.weightA, self.weightB], 
                    [self.proA, self.proB], self.gop, self.scale, self.factor,
                        self.scorer, '1', mode)
            assert alignments[0][-1] < 1
            alignments = _calign.align_pairwise([self.seqA2, self.seqB2],
                    [self.weightA2, self.weightB2], 
                    [self.proA2, self.proB2], self.gop, self.scale, self.factor,
                        self.scorer, '1', mode)
            assert alignments[0][-1] < 1 
            alignments = _calign.align_pairs([[self.seqA, self.seqB]],
                    [[self.weightA, self.weightB]], 
                    [[self.proA, self.proB]], self.gop, self.scale, self.factor,
                        self.scorer, mode, '1', 1)
            assert alignments[0][-1] < 1
            alignments = _calign.align_pairs([[self.seqA2, self.seqB2]],
                    [[self.weightA2, self.weightB2]], 
                    [[self.proA2, self.proB2]], self.gop, self.scale, self.factor,
                        self.scorer, mode, '1', 0)
            assert alignments[0][-1] > 1

        for mode in ['global', 'overlap', 'dialign']:
            profiles = _calign.align_profile(self.profileA, self.profileB,
                    self.weightA, self.weightB, self.proA, self.proB, self.gop,
                    self.scale, self.factor, self.scorer, '1', mode, 0.5)
            assert profiles[-1] == 0.0
            profiles = _calign.align_profile(self.profileA2, self.profileB2,
                    self.weightA2, self.weightB2, self.proA2, self.proB2, self.gop,
                    self.scale, self.factor, self.scorer, '1', mode, 0.5)
            assert len(profiles) == 3


        assert _calign.score_profile(
                ['a', 'a'], ['a', 'a'], self.scorer, 0
                ) == 1

        assert _calign.swap_score_profile(
                ['a', '+'], ['X', 'X'], self.scorer, 0
                ) < 0




    def test_corrdist(self):

        threshold = 0.5
        seqs1 = [[self.seqA, self.seqB]]
        seqs2 = [[self.seqA2, self.seqB2]]
        gops1 = [[self.weightA, self.weightB]]
        gops2 = [[self.weightA2, self.weightB2]]
        pros1 = [[self.proA, self.proB]]
        pros2 = [[self.proA2, self.proB2]]

        for mode in ['global', 'local', 'dialign', 'overlap']:
            
            corr1 = _calign.corrdist(0.5, seqs1, gops1, pros1, self.gop,
                    self.scale,
                    self.factor, self.scorer, mode, '1')
            corr2 = _calign.corrdist(0.5, seqs2, gops2, pros2, self.gop,
                    self.scale,
                    self.factor, self.scorer, mode, '1')
            assert corr1[0]['b', 'b'] == 2
            assert corr2[0]['a', 'a'] == 2 



       
