# *-* coding: utf-8 *-*

"""
Module providing various methods for using Ngram models.

The smoothing methods are implemented to be as compatible as possible with
those offered by NLTK. In fact, both implementation and comments try to follow
Bird at al. as close as possible.
"""

from __future__ import division
import math
import random
from functools import partial

# Try to load the scientific libraries needed for Simple Good-Turing smoothing.
try:
    import numpy as np
except ImportError:
    np = False

try:
    from scipy import linalg, stats
except ImportError:
    linalg, stats = False, False

# Default probability for unobserved samples.
_UNOBS = 1e-10

def _check_probdist_args(freqdist, **kwargs):
    """
    Internal function for validing arguments for smoothing functions.

    Not intended to be called directly by users.
    """

    # Make sure we have a dictionary for the frequency distribution.
    if not isinstance(freqdist, dict):
        raise ValueError("Frequency distribution must be a dictionary.")

    # Get arguments, one by one, and check them; we default to None so we can
    # skip over in case an argument was not provided.
    unobs_prob = kwargs.get('unobs_prob', None)
    if unobs_prob:
        if unobs_prob < 0.0 or unobs_prob > 1.0:
            raise ValueError("Reserved unobserved probability must be in "
                             "range [0.0, 1.0].")

    gamma = kwargs.get('gamma', None)
    if gamma:
        if gamma < 0:
            raise ValueError("Gamma must be a real number.")

    bins = kwargs.get('bins', None)
    if bins:
        if bins < 0:
            raise ValueError("Number of bins must be a real number.")

    default_p0 = kwargs.get('default_p0', None)
    if default_p0:
        if default_p0 < 0 or default_p0 > 1:
            raise ValueError("P0 probability must be in range [0.0, 1.0].")

    p_value = kwargs.get('p_value', None)
    if p_value:
        if p_value <= 0.0 or p_value >= 1.0:
            raise ValueError("p-value must be in range (0.0, 1.0).")


# This kind of work-around to keeping track of which smoothing method was used
# (needed for easier serialization) is not the most elegant or efficient
# (particularly when considering the if/elif structure here employed), but
# should be acceptable for our purposes in order to have a simpler interface
# for users.
# NOTE: Be sure to update the docstring for this function and the one
#       for NgramModel.train() in case of changes to the available smoothing
#       methods.
def smooth_dist(freqdist, method, **kwargs):
    """
    Returns a smoothed log-probability distribution from a named method.

    This method is used to generalize over all implemented smoothing methods,
    especially in terms of serialization. The `method` argument informs which
    smoothing mehtod to use and passes all the arguments to the appropriate
    function.

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the log-probability distribution will be calculated.

    method: str
        The name of the probability smoothing method to use. Either "uniform",
        "random", "mle", "lidstone", "laplace", "ele", "wittenbell",
        "certaintydegree", or "sgt".

    kwargs: additional arguments
        Additional arguments passed to the appropriate smoothing method
        function.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in
        the frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.

    """

    if method == 'uniform':
        sm_func = uniform_dist
    elif method == 'random':
        sm_func = random_dist
    elif method == 'mle':
        sm_func = mle_dist
    elif method == 'lidstone':
        sm_func = lidstone_dist
    elif method == 'laplace':
        sm_func = laplace_dist
    elif method == 'ele':
        sm_func = ele_dist
    elif method == 'wittenbell':
        sm_func = wittenbell_dist
    elif method == 'certaintydegree':
        sm_func = certaintydegree_dist
    elif method == 'sgt':
        sm_func = sgt_dist
    else:
        raise ValueError("Unknown probability smoothing method '%s'." % method)

    return sm_func(freqdist, **kwargs)

def uniform_dist(freqdist, **kwargs):
    """
    Returns a uniform log-probability distribution.

    In a uniform log-probability distribution all samples, no matter the
    observed counts, will have the same log-probability. A mass probability can
    optionally be reserved for unobserved samples.

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the log-probability distribution will be calculated.

    unobs_prob : float
        An optional mass probability to be reserved for unobserved states,
        from 0.0 to 1.0.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in
        the frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Deal with additional arguments.
    unobs_prob = kwargs.get('unobs_prob', _UNOBS)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, unobs_prob=unobs_prob)

    # Calculation couldn't be easier: we just substract the reserved mass
    # probability from 1.0 and divide by the number of samples.
    prob = math.log((1.-unobs_prob) / len(freqdist))

    return {sample:prob for sample in freqdist}, math.log(unobs_prob)


def random_dist(freqdist, **kwargs):
    """
    Returns a random log-probability distribution.

    In a random log-probability distribution all samples, no matter the
    observed counts, will have a random log-probability computed from a set of
    randomly drawn floating point values. A mass probability can optionally be
    reserved for unobserved samples.

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    unobs_prob : float
        An optional mass probability to be reserved for unobserved states,
        from 0.0 to 1.0.

    seed : any hasheable value
        An optional seed for the random number generator, defaulting to None.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in
        the frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Deal with additional arguments.
    unobs_prob = kwargs.get('unobs_prob', _UNOBS)
    seed = kwargs.get('seed', None)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, unobs_prob=unobs_prob)

    # Set the seed, generate a random number for each sample and the sum of
    # fake observations, and build a probability distribution. We sort the keys
    # in the frequency distribution to guarantee that a set of samples will
    # always have the same random distribution (especially if a seed is
    # provided) -- as `freqdist` is a dictionary, in most Python
    # implementations the order of the keys is not guaranteed.
    random.seed(seed)
    fake_count = {sample:random.random() for sample in sorted(freqdist)}
    fake_sum = sum(fake_count.values())

    probdist = {sample:math.log((count / fake_sum) * (1.-unobs_prob))
                for sample, count in fake_count.items()}

    return probdist, math.log(unobs_prob)


def mle_dist(freqdist, **kwargs):
    """
    Returns a Maximum-Likelihood Estimation log-probability distribution.

    In an MLE log-probability distribution the probability of each sample is
    approximated as the frequency of the same sample in the frequency
    distribution of observed samples. It is the distribution people intuitively
    adopt when thinking of probability distributions. A mass probability can
    optionally be reserved for unobserved samples.

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    unobs_prob : float
        An optional mass probability to be reserved for unobserved states,
        from 0.0 to 1.0.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Deal with additional arguments.
    unobs_prob = kwargs.get('unobs_prob', _UNOBS)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, unobs_prob=unobs_prob)

    # Run the estimator by simply collecting the sum of values and dividing the
    # counts of each sample by such value.
    value_sum = sum(freqdist.values())
    probdist = {sample:math.log((count/value_sum) * (1.-unobs_prob))
                for sample, count in freqdist.items()}

    return probdist, math.log(unobs_prob)


def lidstone_dist(freqdist, **kwargs):
    """
    Returns a Lidstone estimate log-probability distribution.

    In a Lidstone estimate log-probability the frequency distribution of
    observed samples is used to estimate the probability distribution of the
    experiment that generated such observation, following a parameter given by
    a real number *gamma* typycally randing from 0.0 to 1.0. The Lidstone
    estimate approximates the probability of a sample with count *c* from an
    experiment with *N* outcomes and *B* bins as *(c+gamma)/(N+B*gamma)*. This
    is equivalent to adding *gamma* to the count of each bin and taking the
    Maximum-Likelihood estimate of the resulting frequency distribution, with
    the corrected space of observation; the probability for an unobserved
    sample is given by frequency of a sample with gamma observations.

    Also called "additive smoothing", this estimation method is frequently
    used with a *gamma* of 1.0 (the so-called "Laplace smoothing") or of 0.5
    (the so-called "Expected likelihood estimate", or ELE).

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    gamma : float
        A real number used to parameterize the estimate.

    bins: int
        The optional number of sample bins that can be generated by the
        experiment that is described by the probability distribution. If not
        specified, it will default to the number of samples in the frequency
        distribution.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Deal with additional arguments.
    gamma = kwargs.get('gamma', None)
    bins = kwargs.get('bins', None)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, gamma=gamma, bins=bins)

    # Obtain the parameters for probability calculation.
    N = sum(freqdist.values())
    if not bins:
        B = len(freqdist)
    else:
        B = bins

    probdist = {sample:math.log((count + gamma) / (N + B * gamma))
                for sample, count in freqdist.items()}
    prob_unk = math.log(gamma / (N + B * gamma))

    return probdist, prob_unk


laplace_dist = partial(lidstone_dist, gamma=1)
laplace_dist.__doc__ = """
    Returns a Laplace estimate log-probability distribution.
    
    In a Laplace estimate log-probability the frequency distribution of
    observed samples is used to estimate the probability distribution of the
    experiment that generated such observation, following a parameter given by
    a real number *gamma* set by definition to 1. As such, it is a
    generalization of the Lidstone estimate.
    
    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    bins: int
        The optional number of sample bins that can be generated by the
        experiment that is described by the probability distribution. If
        not specified, it will default to the number of samples in
        the frequency distribution.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """


ele_dist = partial(lidstone_dist, gamma=0.5)
ele_dist.__doc__ = """
    Returns an Expected-Likelihood estimate log-probability distribution.
    
    In an Expected-Likelihood estimate log-probability the frequency
    distribution of observed samples is used to estimate the probability
    distribution of the experiment that generated such observation, following a
    parameter given by a real number *gamma* set by definition to 0.5. As such,
    it is a generalization of the Lidstone estimate.
    
    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    bins: int
        The optional number of sample bins that can be generated by the
        experiment that is described by the probability distribution. If not
        specified, it will default to the number of samples in the frequency
        distribution.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """


def wittenbell_dist(freqdist, **kwargs):
    """
    Returns a Witten-Bell estimate log-probability distribution.

    In a Witten-Bell estimate log-probability a uniform probability mass is
    allocated to yet unobserved samples by using the number of samples that
    have only been observed once. The probability mass reserved for unobserved
    samples is equal to *T / (N +T)*, where *T* is the number of observed
    samples and *N* the number of total observations. This equates to the
    Maximum-Likelihood Estimate of a new type of sample occurring. The
    remaining probability mass is discounted such that all probability
    estimates sum to one, yielding:

        - *p = T / Z (N + T)*, if count == 0
        - *p = c / (N + T)*, otherwise

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    bins: int
        The optional number of sample bins that can be generated by the
        experiment that is described by the probability distribution. If not
        specified, it will default to the number of samples in the frequency
        distribution.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Deal with additional arguments.
    bins = kwargs.get('bins', None)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, bins=bins)

    # Obtain the parameters for probability calculation; we are replacing `B`
    # by `T` as a notation, here, to make it clear that it is not necessarily
    # the same `B` value computed in other probability distributions.
    N = sum(freqdist.values())
    T = len(freqdist)
    if not bins:
        Z = 1.0
    else:
        if bins == T:
            Z = 1.0
        else:
            Z = bins - T

    # Compute the log-probabilities.
    probdist = {sample:math.log(count / (N + T))
                for sample, count in freqdist.items()}

    # The probability for unobserved samples depends on N.
    if N == 0:
        prob_unk = math.log(1.0 / Z)
    else:
        prob_unk = math.log(T / (Z * (N + T)))

    return probdist, prob_unk


def certaintydegree_dist(freqdist, **kwargs):
    """
    Returns a log-probability distribution based on the degree of certainty.

    In this distribution a mass probability is reserved for unobserved samples
    from a computation of the degree of certainty that the are no unobserved
    samples.

    Under development and test by Tiago Tresoldi, this is an experimental
    probability distribution that should not be used as the sole or main
    distribution for the time being.

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    bins: int
        The optional number of sample bins that can be generated by the
        experiment that is described by the probability distribution. If not
        specified, it will default to the number of samples in the frequency
        distribution.

    unobs_prob : float
        An optional mass probability to be reserved for unobserved states,
        from 0.0 to 1.0.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Deal with additional arguments.
    bins = kwargs.get('bins', None)
    unobs_prob = kwargs.get('unobs_prob', _UNOBS)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, bins=bins)

    # Obtain the parameters for probability calculation.
    N = sum(freqdist.values())
    B = len(freqdist)
    Z = bins or B

    # Calculate the mass of probability space to reserve and use this value to
    # correct the Maximum-Likelihood Estimate for each sample.
    # NOTE: For very large values of N, this will underflow because we
    #       effectively have a large confidence of having observed all the
    #       samples that matter; this is a problem when taking the
    #       log-probability, as we'll ultimately raise a math domain error by
    #       asking for the logarithm of what is machine-represented as zero;
    #       for this reason, we take as the probability space the minimum value
    #       between 1.0 discounted the calculated mass and 1.0 discounted the
    #       minimum mass probability reserved.
    prob_space = min(1. - (B/(Z+1))**N, 1. - unobs_prob)
    probdist = {sample:math.log((count/N) * prob_space)
                for sample, count in freqdist.items()}

    # Calculate the unobserved probability from the probability space.
    prob_unk = math.log(-(prob_space - 1.))

    return probdist, prob_unk

def sgt_dist(freqdist, **kwargs):
    """
    Returns a Simple Good-Turing log-probability distribution.

    The returned log-probability distribution is based on the Good-Turing
    frequency estimation, as first developed by Alan Turing and I. J. Good and
    implemented in a more easily computable way by Gale and Sampson's
    (1995/2001 reprint) in the so-called "Simple Good-Turing".

    This implementation is based mostly in the one by "maxbane" (2011)
    (https://github.com/maxbane/simplegoodturing/blob/master/sgt.py), as well
    as in the original one in C by Geoffrey Sampson (1995; 2000; 2005; 2008)
    (https://www.grsampson.net/Resources.html), and in the one by
    Loper, Bird et al. (2001-2018, NLTK Project)
    (http://www.nltk.org/_modules/nltk/probability.html). Please note that
    due to minor differences in implementation intended to guarantee non-zero
    probabilities even in cases of expected underflow, as well as our
    relience on scipy's libraries for speed and our way of handling
    probabilities that are not computable when the assumptions of SGT are
    not met, most results will not exactly match those of the 'gold standard'
    of Gale and Sampson, even though the differences are never expected to
    be significative and are equally distributed across the samples.

    Parameters
    ----------

    freqdist : dict
        Frequency distribution of samples (keys) and counts (values) from
        which the probability distribution will be calculated.

    p_value : float
        The p-value for calculating the confidence interval of the empirical
        Turing estimate, which guides the decision of using either the Turing
        estimate "x" or the loglinear smoothed "y". Defaults to 0.05, as per
        the reference implementation by Sampson, but consider that the authors,
        both in their paper and in the code following suggestions credited to
        private communication with Fan Yang, consider using a value of 0.1.

    allow_fail : bool
        A logic value informing if the function is allowed to fail, throwing
        RuntimeWarning exceptions, if the essential assumptions on the
        frequency distribution are not met, i.e., if the slope of the loglinear
        regression is > -1.0 or if an unobserved count is reached before we are
        able to cross the smoothing threshold. If set to False, the estimation
        might result in an unreliable probability distribution; defaults to
        True.

    default_p0 : float
        An optional value indicating the probability for unobserved samples
        ("p0") in cases where no samples with a single count are observed; if
        this value is not specified, "p0" will default to a Laplace estimation
        for the current frequency distribution. Please note that this is
        intended change from the reference implementation by Gale and Sampson.

    Returns
    -------

    state_prob: dict
        A dictionary of sample to log-probabilities for all the samples in the
        frequency distribution.

    unobserved_prob: float
        The log-probability for samples not found in the frequency
        distribution.
    """

    # Make sure the scientific libraries have been loaded, raising an
    # ImportError if not
    if not np:
        raise ImportError('The package `numpy` is needed by SGT.')
    if not linalg or not stats:
        raise ImportError('The package `scipy` is needed by SGT.')

    # Deal with additional arguments.
    default_p0 = kwargs.get('default_p0', None)
    p_value = kwargs.get('p_value', 0.05)
    allow_fail = kwargs.get('allow_fail', True)

    # Perform basic argument checking.
    _check_probdist_args(freqdist, default_p0=default_p0, p_value=p_value)

    # Calculate the confidence level from the p_value.
    confidence_level = stats.norm.ppf(1. - (p_value/2.0))

    # Remove all samples with `count` equal to zero.
    freqdist = {sample:count
                for sample, count in freqdist.items() if count > 0}

    # Prepare vectors for frequencies (`r` in G&S) and frequencies of
    # frequencies (`Nr` in G&S). freqdist.values() is cast to a tuple because
    # we can't consume the iterable a single time. `freqs_keys` is sorted to
    # make vector computations faster later on (so we query lists and not
    # dictionaries).
    freqs = tuple(freqdist.values())
    freqs_keys = sorted(set(freqs)) # r -> n (G&S)
    freqs_of_freqs = {c:freqs.count(c) for c in freqs_keys}

    # The papers and the implementations are not clear on how to calculate the
    # probability of unobserved states in case of missing single-count samples
    # (unless we just fail, of course); Gale and Sampson's C implementation
    # defaults to 0.0, which is not acceptable for our purposes. The solution
    # here offered is to either use an user-provided probability (but in this
    # case we are not necessarily defaulting to _UNOBS, and, in fact, the
    # function argument name is `default_p0` and not `unobs_prob`) or default
    # to a Lidstone smoothing with a gamma of 1.0 (i.e., using Laplace
    # smoothing constant).
    # TODO: Investigate and discuss other possible solutions, including
    #       user-defined `gamma`, `bins`, and/or `N`.
    if 1 in freqs_keys:
        p0 = freqs_of_freqs[1] / sum(freqs)
    else:
        p0 = default_p0 or (1. / (sum(freqs)+1))

    # Compute Sampson's Z: for each count `j`, we set Z[j] to the linear
    # interpolation of {i, j, k}, where `i` is the greatest observed count less
    # than `j`, and `k` the smallest observed count greater than `j`.
    I = [0] + freqs_keys[:-1]
    K = freqs_keys[1:] + [2*freqs_keys[-1] - I[-1]]
    Z = {j:2*freqs_of_freqs[j] / (k-i)
         for i, j, k in zip(I, freqs_keys, K)}

    # Compute a loglinear regression of Z[r] over r. We cast keys and values to
    # a list for the computation with `linalg.lstsq`.
    z_keys = list(Z.keys())
    z_values = list(Z.values())
    slope, intercept = \
        linalg.lstsq(np.c_[np.log(z_keys), (1,)*len(z_keys)],
                     np.log(z_values))[0]
    #print ('Regression: log(z) = %f*log(r) + %f' % (slope, intercept))
    if slope > -1.0 and allow_fail:
        raise RuntimeWarning("In SGT, linear regression slope is > -1.0.")

    # Aapply Gale and Sampson's "simple" loglinear smoothing method.
    r_smoothed = {}
    use_y = False
    for r in freqs_keys:
        # `y` is the loglinear smoothing.
        y = float(r+1) * \
            np.exp(slope*np.log(r+1) + intercept) / \
            np.exp(slope*np.log(r) + intercept)

        # If we've already started using `y` as the estimate for `r`, then
        # continue doing so; also start doing so if no samples were observed
        # with count equal to `r+1` (following comments and variable names in
        # both Sampson's C implementation and in NLTK, we check at which
        # point we should `switch`)
        if r+1 not in freqs_of_freqs:
            if not use_y:
                # An unobserved count was reached before we were able to cross
                # the smoothing threshold; this means that assumptions were
                # not met and the results will likely be off.
                if allow_fail:
                    raise RuntimeWarning("In SGT, unobserved count before smoothing threshold.")

            use_y = True

        # If we are using `y`, just copy its value to `r_smoothed`, otherwise
        # perform the actual calculation.
        if use_y:
            r_smoothed[r] = y
        else:
            # `estim` is the empirical Turing estimate for `r` (equivalent to
            # `x` in G&S)
            estim = (float(r+1) * freqs_of_freqs[r+1]) / freqs_of_freqs[r]

            Nr = float(freqs_of_freqs[r])
            Nr1 = float(freqs_of_freqs[r+1])

            # `width` is the width of the confidence interval of the empirical
            # Turing estimate (for which Sampson uses 95% but suggests 90%),
            # when assuming independence.
            width = confidence_level * \
                    np.sqrt(float(r+1)**2 * (Nr1 / Nr**2) * (1. + (Nr1 / Nr)))

            # If the difference between `x` and `y` is more than `t`, then the
            # empirical Turing estimate `x` tends to be more accurate.
            # Otherwise, use the loglinear smoothed value `y`.
            if abs(estim - y) > width:
                r_smoothed[r] = estim
            else:
                use_y = True
                r_smoothed[r] = y

    # (Re)normalize and return the resulting smoothed probabilities, less the
    # estimated probability mass of unseen species; please note that we might
    # be unable to calculate some probabilities if the function was not allowed
    # to fail, mostly due to math domain errors. We default to `p0` in all such
    # cases.
    smooth_sum = sum([freqs_of_freqs[r] * r_smooth
                      for r, r_smooth in r_smoothed.items()])

    # Build the probability distribution for the observed samples and for
    # unobserved ones.
    prob_unk = math.log(p0)
    probdist = {}
    for sample, count in freqdist.items():
        prob = (1.0 - p0) * (r_smoothed[count] / smooth_sum)
        if prob == 0.0:
            probdist[sample] = math.log(p0)
        else:
            probdist[sample] = math.log(prob)

    return probdist, prob_unk
