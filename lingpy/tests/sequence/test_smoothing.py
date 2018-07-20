# *-* coding: utf-8 *-*

import itertools
import random
from collections import Counter
from nose.tools import assert_raises
from unittest import TestCase

from lingpy.sequence.smoothing import smooth_dist

"""
Tests for the smoothing module.
"""

class Tests(TestCase):
    def setUp(self):
        # Set up a bunch of fake distribution to test the methods.
        # We don't use real data as this is only intended to test
        # the programming side, and we don't want to distribute a
        # datafile only for this testing purposes. We setup a simple
        # distribution with character from A to E, and a far more
        # complex one (which includes almost all printable characters
        # as single states, plus a combination of letters and characters,
        # all with a randomly generated frequency)
        self.observ1 = Counter([char for char in 'ABBCCCDDDDEEEE'])

        samples = [char for char in
            '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$%&()*+,-./:;?@[\\]^_`{|}~']
        two_char_samples = [[char1+char2 for char2 in '0123456789'] for char1 in 'ABCDEFGHIJK']
        samples += itertools.chain.from_iterable(two_char_samples)
        
        random.seed(1305)
        self.observ2 = {
            sample:(random.randint(1, 1000)**random.randint(1, 3))+random.randint(1, 100)
            for sample in samples}
        self.observ2['x'] = 100
        self.observ2['y'] = 100
        self.observ2['z'] = 1000

    def test_uniform_dist(self):
        """
        Test for the uniform distribution.
        """

        # Easiest distribution to test: everything must be equal, unobserved
        # must be less than any observed one.
        seen, unseen = smooth_dist(self.observ1, "uniform")
        assert seen['A'] == seen['E']
        assert seen['A'] > unseen

        seen, unseen = smooth_dist(self.observ2, "uniform")
        assert seen['0'] == seen['~']
        assert seen['x'] == seen['z']
        assert seen['0'] > unseen


    def test_random_dist(self):
        """
        Test for the random distribution.
        """

        # Here, we can essentially only test if there are no execution bugs,
        # besides the unobserved probability being less than any other. Nonetheless,
        # we set a seed for comparing the results across platforms.
        seen, unseen = smooth_dist(self.observ1, "random", seed=1305)
        assert seen['A'] < seen['E']
        assert seen['B'] < seen['D']
        assert seen['A'] < seen['C']
        assert seen['A'] > unseen

        seen, unseen = smooth_dist(self.observ2, "random", seed=1305)
        assert seen['B1'] < seen['4']
        assert seen['F5'] < seen['C8']
        assert seen['0'] > unseen

    def test_mle_dist(self):
        """
        Test for the Maximum-Likelihood Estimation distribution.
        """

        # This is easy to test as the results of the MLE are the ones we intuitively
        # expect/compute.
        seen, unseen = smooth_dist(self.observ1, "mle")
        assert seen['A'] < seen['B']
        assert seen['D'] == seen['E']
        assert seen['A'] > unseen

        seen, unseen = smooth_dist(self.observ2, "mle")
        assert seen['x'] == seen['y']
        assert seen['x'] < seen['z']
        assert seen['0'] > unseen

    def test_laplace_dist(self):
        """
        Test for the Laplace distribution.
        """

        seen, unseen = smooth_dist(self.observ1, "laplace")
        assert seen['A'] < seen['B']
        assert seen['D'] == seen['E']
        assert seen['A'] > unseen

        seen, unseen = smooth_dist(self.observ2, "laplace")
        assert seen['x'] == seen['y']
        assert seen['x'] < seen['z']
        assert seen['0'] > unseen

    def test_ele_dist(self):
        """
        Test for the Expected-Likelihood estimation distribution.
        """

        seen, unseen = smooth_dist(self.observ1, "ele")
        assert seen['A'] < seen['B']
        assert seen['D'] == seen['E']
        assert seen['A'] > unseen

        seen, unseen = smooth_dist(self.observ2, "ele")
        assert seen['x'] == seen['y']
        assert seen['x'] < seen['z']
        assert seen['0'] > unseen

    def test_wittenbell_dist(self):
        """
        Test for the Witten-Bell distribution.
        """

        seen10, unseen10 = smooth_dist(self.observ1, "wittenbell", bins=10)
        seen99, unseen99 = smooth_dist(self.observ1, "wittenbell", bins=99)
        assert seen10['A'] < seen10['B']
        assert seen10['D'] == seen10['E']
        assert seen10['A'] == unseen10
        assert seen10['A'] == seen99['A']
        assert unseen99 < unseen10

        seen, unseen = smooth_dist(self.observ2, "wittenbell")
        assert seen['x'] == seen['y']
        assert seen['x'] < seen['z']
        assert seen['0'] > unseen

    def test_certaintydegree_dist(self):
        """
        Test for the Degree of Certainty distribution.
        """

        seen, unseen = smooth_dist(self.observ1, 'certaintydegree')
        seen99, unseen99 = smooth_dist(self.observ1, 'certaintydegree', bins=99)
        assert seen['A'] < seen['B']
        assert seen['D'] == seen['E']
        assert seen['A'] < unseen
        assert seen['B'] > unseen
        assert seen['A'] < seen99['A']
        assert unseen > unseen99

        seen, unseen = smooth_dist(self.observ2, "certaintydegree")
        assert seen['x'] == seen['y']
        assert seen['x'] < seen['z']
        assert seen['0'] > unseen

    def test_sgt_dist(self):
        """
        Test for the Simple Good-Turing distribution.
        """

        # Only run the test if the numpy and scipy libraries (imported
        # by the function) are installed; an ImportError will not be a
        # test failure as we decided that they should not be dependencies.

        try:
            # The first distribution does not have enough data for SGT,
            # assert that an assertion is raised.
            assert_raises(RuntimeWarning, smooth_dist, self.observ1, "sgt")

            # The second distribution also does not have enough data for
            # a confident results, but it does have enough to stress-test the
            # method, so we are not going to allow it to fail (i.e., raise
            # exceptions)
            seen_p05, unseen_p05 = smooth_dist(self.observ2, "sgt", allow_fail=False)
            assert seen_p05['x'] == seen_p05['y']
            assert seen_p05['x'] < seen_p05['z']
            assert seen_p05['0'] > unseen_p05
        except ImportError:
            pass
