# *-* coding: utf-8 *-*
"""
This modules provides methods for generating and collecting ngrams.

The methods allow to collect different kind of subsequences, such as standard
ngrams (preceding context), skip ngrams with both single or multiple gap
openings (both preceding and following context), and positional ngrams (both
preceding and following context).
"""

from collections import defaultdict, Counter
from itertools import chain, combinations, product
from functools import partial

# Global padding symbol, shared across all functions/class-methods.
PAD_SYMBOL = '$$$'
ELM_SYMBOL = '###'

def _seq_as_tuple(sequence):
    """
    Internal function for automatically converting a string sequence to a
    tuple, if needed.

    Parameters
    ----------
    sequence: list or str
        The sequence that shall be converted into an iterable.

    Returns
    -------
    out: tuple
        A tuple of the sequence.
    """

    # We first check for datatype and then for a space, as the first test is
    # faster (and evaluation is lazy).
    if isinstance(sequence, str) and ' ' in sequence:
        return tuple(sequence.split(' '))

    return tuple(sequence)


def _random_choices(population, weights=None, cum_weights=None, k=1):
    """
    Return a `k` sized list of elements chosen from `population` with
    replacement and according to a list of weights.

    If a `weights` sequence is specified, selections are made according to the
    relative weights. Alternatively, if a `cum_weights` sequence is given, the
    selections are made according to the cumulative weights. For example, the
    relative weights `[10, 5, 30, 5]` are equivalent to the cumulative weights
    `[10, 15, 45, 50]`. Internally, the relative weights are converted to the
    cumulative weights before making selections, so supplying the cumulative
    weights saves work.

    This function is compatible with the random.choices() function available
    in Python's standard library from version 3.6 on. It can be replaced by
    the standard implementation once the version requirement is updated.

    Parameters
    ----------
    population: list
        A list of elements from which the element(s) will be drawn.

    weights: list
        A list of any numeric type with the relative weight of each element.
        Either `weights` or `cum_weights` must be provided.

    cum_weights: list
        A list of any numeric type with the accumulated weight of each element.
        Either `weights` or `cum_weights` must be provided.

    k: int
        The number of elements to be drawn, with replacement.

    Returns
    -------
    sample: list
        A list of elements randomly drawn according to the specified weights.
    """

    # Assert that (1) the population is not empty, (2) only one type of
    # weight information is provided.
    assert len(population) > 0, "Population must not be empty."
    assert not all((weights, cum_weights)), \
        "Either only weights or only cumulative weights must be provided."

    # If cumulative weights were not provided, build them from `weights`.
    if not cum_weights:
        cum_weights = list(itertools.accumulate(weights))

    # Assert that the lengths of population and cumulative weights match.
    assert len(population)==len(cum_weights), \
        "Population and weight lengths do not match."

    # Get a random number and see in which bin it falls. We need to use this
    # logic which is a little more complex than something with randint()
    # in order to allow for floating-point weights.
    rnd = [random.uniform(0, cum_weights[-1]) for r in range(k)]
    less_than = [[cw < r for cw in cum_weights] for r in rnd]

    return [population[lt.index(False)] for lt in less_than]


class NgramModel():
    """
    Class for operation upon sequences using ngrams models.

    This class allows different operations upon sequences after training ngram
    models, such as sequence relative likelihood computation (both per state
    and overall), random sequence generation, computation of a model entropy
    and of cross-entropy/perplexity of a sequence, etc. As model training is
    computationally and time consuming for large datasets, trained models can
    be saved and loaded ("serialized") from disk.
    """

    def __init__(self, pre_order=None, post_order=None, pad_symbol=PAD_SYMBOL,
                 sequences=[], model_file=None):
        """
        Add docstring
        """

        # Initialize internal parameters from the arguments. Users are
        # allowed to pass either a list of sequences or a filename with a
        # serialized model; in case both are informed, the model is loaded
        # first (including the ngram counting), the smoothed probabilities are
        # cleaned, and the model is update with the new sequences.

        # Store the ngram collection parameters. While the user can pass
        # `pre_order` and `post_order` as integers (indicating the maximum
        # order) or lists (allowing to specify discontinuous distirbutions),
        # just as in the functions for collecting ngrams, we are already
        # casting the range iterators to lists here, so we can quickly query
        # the values; this is necessary, for example, for random sequence
        # generation when we need to find the largest order. This behaviour
        # should be entirely transparent to users.
        self._padsymbol = pad_symbol

        if isinstance(pre_order, int):
            self._pre = list(range(pre_order + 1))
        else:
            self._pre = pre_order

        if isinstance(post_order, int):
            self._post = list(range(post_order + 1))
        else:
            self._post = post_order

        # Initialize internal variables for holding the model.
        self._ngrams = defaultdict(lambda: Counter())
        self._ngram_space = Counter()
        self._seqlens = Counter()

        # The number of bins that was used/will be used in training.
        self._bins = None

        # Whether the model was trained.
        self._trained = False

        # If a model file was provided but no new sequences, just load the
        # model file (which already includes the trained model). If only
        # sequences were provided but no model file, collect the ngrams from
        # the sequences and don't perform any training yet. If both a model
        # file was provided and new sequences, load the model file, clean
        # the training (as we have new sentences), and collected the new
        # ngrams (no automatic training is performed).
        if model_file:
            self.load_model(model_file)

        # Call add_sequences(); a check of whether sequences were provided is
        # performed in the function itself.
        self.add_sequences(sequences)

    def add_sequences(self, sequences):
        """
        Adds sequences to a model, collecting their ngrams.

        This method does not return any value, but cleans the internal
        matrix probability, if one was previously computed, and automatically
        updates the ngram counters. The actual training, with the computation
        of smoothed log-probabilities, is not performed automatically, and
        must be requested by the user by calling the .train() method.

        Parameters
        ----------
        sequences: list
            A list of sequences to be added to the model.
       """

        if sequences:
            # Either initialize (if no model file was provided) or clear (if
            # a model file was provided) the variables for smoothed
            # probabilities.
            self._p = {}
            self._p0 = {}
            self._l = {}
            self._l0 = {}
            self._trained = False

            for sequence in sequences:
                # Collect all positional ngrams, using the ngram tuple as a key
                # and the state as value (which is appended to self._ngrams()).
                # The positional information (ngram[2]) is actually discarded
                # in this stage.
                [self._ngrams[ngram[0]].update(ngram[1]) for ngram in
                 get_all_posngrams(sequence, self._pre, self._post, self._padsymbol)]

            # Collect sequence lengths.
            self._seqlens.update([len(sequence) for sequence in sequences])

    def train(self, method='laplace', normalize=False, bins=None):
        """
        Train a model after ngrams have been collected.

        This method does not return any value, but sets the internal variables
        with smoothed probabilities (such as `self._p` and `self._p0`) and
        internally marks the model as having been trained.

        Parameters
        ----------
        method: str
            The name of the smoothing method to be used, as used by
            `smooth_dist()`. Either "uniform", "random", "mle", "lidstone",
            "laplace", "ele", "wittenbell", "certaintydegree", or "sgt".
            Defaults to "laplace".

        normalize: boolean
            Whether to normalize the log-probabilities for each ngram in the
            model after smoothing, i.e., guarantee that the probabilities (with
            the probability for unobserved transitions counted a single time)
            sum to 1.0. This is computationally expansive, and should be only
            used if the model is intended for later serialization. While
            experiments with real data demonstrated that this normalization
            does not improve the results or performance of the methods, the
            computational cost of normalizing the probabilities might be
            justified if descriptive statistics on the model, like samples
            from the matrix of transition probabilities or the
            entropy/perplexity of a sequence, are needed (such as for
            publication), as they will be more in line with what is generally
            expected and will facilitate the comparison of different models.

        bins: int
            The number of bins to be assumed when smoothing, for the smoothing
            methods that use this information. Defaults to the number of
            unique states observed, as gathered from the count of ngrams with
            no context.
       """

        # No need to initialize/clean `self._p` and `self._p0` (as well as the
        # corresponding variables for length probabilities) as they are
        # cleaned by `self.add_sequences()` everytime a sequence is added.

        # If the number of bins was not informed, use the number of transition
        # states from the zero-context.
        # TODO: assert that we have at least the number of observed transitions?
        # TODO: what to do if there is no no-context ngram? Should be
        # mandatory?
        self._bins = bins or len(self._ngrams[('###',)])

        # Perform the probability smoothing.
        # TODO: allow correction in which all observed probs smoothed to less
        # than non-observed are equal to non-observed at least -- should this
        # be moved to smooth_dist?
        for context, counter in self._ngrams.items():
            self._p[context], self._p0[context] = \
                smooth_dist(counter, method=method, bins=self._bins)

        # Normalize, if so requested. See comments in the docstring for more
        # information.
        if normalize:
            for context, probs in self._p.items():
                # We remap the log-probabilities into probabilities in the
                # temporary variable `_prob`, sum them in order to obtain the
                # used probability space (which included a single occurence
                # for the unobserved probability), and finally reset
                # `self._p` and `self._p0` as normalized log-probabilities.
                _prob = {state:math.exp(prob) for state, prob in probs.items()}
                _prob0 = math.exp(self._p0[context])
                _prob_sum = sum(_prob.values()) + _prob0

                self._p[context] = {state:math.log(prob/prob_sum)
                                    for state, prob in _prob.items()}
                self._p0[context] = math.log(_prob0/prob_sum)

        # Compute the log-probabilities for lengths. This is easy as we just
        # assume that the count/probability for non-observed lengths is equal
        # to the count/probability of the less observed length (the value is
        # added directly to `length_obs`). This is similar to ML estimation.
        length_obs = sum(self._seqlens.values()) + min(self._seqlens.values())
        self._lp = {length:math.log(count/length_obs)
                    for length, count in self._seqlens.items()}
        self._lp0 = math.log(min(self._seqlens.values()) / length_obs)

        # Collect the ngram space keys and values for random sequence
        # generation.
        # TODO: do it at the same time that the other collections
        for context, counter in self._ngrams.items():
            for key, value in counter.items():
                key = tuple(s if s != '###' else key for s in context)
                self._ngram_space[key] += value

        # Internally inform that the model was trained.
        self._trained = True

    # TODO: correct not found by no context probability, optionally
    def state_score(self, sequence):
        """
        Returns the relative likelihood for each state in a sequence.

        Please note that this does not perform correction due to sequence
        length. The model must have been trained before using this function.

        Parameters
        ----------
        sequence: list
            A list of states to be scored.

        Returns
        -------
        prob: list
            A list of floats, of the same length of `sequence`, with the
            individual log-probability for each state.
        """

        # Assert the model was trained.
        assert self._trained==True, "Ngram Model was not trained."

        # Pre-allocate the list holding the probability (i.e., the relative
        # likelihood) for each state in `sequence`.
        s_prob = [0.0] * len(sequence)

        # We collect all positional ngrams in `sequence`, using the same
        # paramenters for the model ngram collection, and compute the state
        # probability in each case, appending/adding the results to the
        # correct element in `s_prob`.

        for ngram, state, idx in \
            get_all_posngrams(sequence, self._pre, self._post, self._padsymbol):
            # If the ngram (the "context") is found in `self._p` (i.e., it was
            # observed in training), we just need to append to `_p` the
            # probability of its state (or of the transition to unonserved
            # states, in case the ngram was observed in training but not the
            # state). If the ngram was not observed, we need a different
            # backoff solution; here we rely in the chain rule and just sum
            # the log-probabilities of each individual state in the ngram
            # (including the state being observed), taking care of
            # unobserved states.
            if ngram in self._p:
                if state in self._p[ngram]:
                    _p = self._p[ngram][state]
                else:
                    # TODO: correction?
                    _p = self._p0[ngram]
            else:
                # Make a copy of the sequence replacing the symbol for the
                # current state by the observed state; then, compute and sum
                # the individual log-probabilities.
                _seq = [state if seq_state=='###' else seq_state
                        for seq_state in ngram]
                _p = sum([self._p[('###',)][seq_state]
                          if seq_state in self._p[('###',)]
                          else self._p0[('###',)]
                          for seq_state in _seq])

            # Update the log-probability of the correct state.
            s_prob[idx] += _p

        return s_prob

    def score(self, sequence, use_length=True):
        """
        Returns the relative likelihood of a sequence.

        The model must have been trained before using this function.

        Parameters
        ----------
        sequence: list
            A list of states to be scored.

        use_length: bool
            Whether to correct the sequence relative likelihood by using
            length probability. Defaults to True.

        Returns
        -------
        prob: list
            A list of floats, of the same length of `sequence`, with the
            individual log-probability for each state.
        """

        # Assert the model was trained.
        assert self._trained==True, "Ngram Model was not trained."

        # Get the sum of individual log-probabilities, correct them with the
        # sequence length probability if requested and return.
        _prob = sum(self.state_score(sequence))
        if use_length:
            if len(sequence) in self._lp:
                _prob += self._lp[len(sequence)]
            else:
                _prob += self._lp0

        return _prob

    def model_entropy(self):
        # collect the P x log(P) for all contexts; probabilities are already
        # stored as ln in self._p, so we need to get back to the probability
        # itself, which is computationally quite expensive
        lentropy = []
        for context in self._p:
            # for manipulation, add _p0 at this stage
            _probs = [math.exp(p) for p in self._p[context].values()]
            _probs.append(math.exp(self._p0[context]))
            lentropy += [-p*math.log2(p) for p in _probs]

        return sum(lentropy)

    def entropy(self, sequence, base=2.0):
        """
        Calculates the cross-entropy of a sequence.

        Parameters
        ----------
        sequence: list
            The sequence whose cross-entropy will be calculated.

        base: float
            The logarithmic base for the cross-entropy calculation. Defaults to
            2.0, following the standard approach set by Shannon that allows
            to consider entropy in terms of bits needed for unique
            representation.

        Returns
        -------
        ce: float
            The cross-entropy calculated for the sequence, a real number.
        """
        return -(self.score(sequence)/math.log(base)) / len(sequence)

    def perplexity(self, sequence):
        """
        Calculates the perplexity of a model.

        As per definition, this is simply 2.0 to the cross-entropy of the
        given sequence on logarithmic base of 2.0.

        Parameters
        ----------
        sequence: list
            The sequence whose perplexity should be calculated.

        Returns
        -------
        perplexity: float
            The calculated perplexity for the sequence.
        """
        return 2.0 ** self.entropy(sequence)

    def percentile_rank(self, sequence):
        raise ValueError("Not implemented yet.")

    # we must use pickle and not json as json's library does not allow for
    # tuple lists, unless we map everythin back and forth from strings
    # while we could iterate over everything, it is better to do it variable
    # by variable, so we take care of special cases and make it easier for
    # future improvements
    def save_model(self, filename):
        # We cannot pickle some data structures, particularly those which use
        # a lambda. For those case, we need to cast back and forth when
        # saving and loading.
        model_data = {
            '_padsymbol' : self._padsymbol,
            '_pre' : self._pre,
            '_post' : self._post,
            '_ngrams' : dict(self._ngrams),
            '_ngram_space' : self._ngram_space,
            '_seqlens' : self._seqlens,
            '_bins' : self._bins,
            '_trained' : self._trained,
            '_p' : self._p,
            '_p0' : self._p0,
            '_l' : self._l,
            '_l0' : self._l0,
        }

        with open(filename, 'wb') as handler:
            pickle.dump(model_data, handler, pickle.HIGHEST_PROTOCOL)

    def load_model(self, filename):
        with open(filename, 'rb') as handler:
            model_data = pickle.load(handler)

        self._padsymbol = model_data['_padsymbol']
        self._pre = model_data['_pre']
        self._post = model_data['_post']
        self._ngrams = defaultdict(lambda: Counter(), model_data['_ngrams'])
        self._ngram_space = model_data['_ngram_space']
        self._seqlens = model_data['_seqlens']
        self._bins = model_data['_bins']
        self._trained = model_data['_trained']
        self._p = model_data['_p']
        self._p0 = model_data['_p0']
        self._l = model_data['_l']
        self._l0 = model_data['_l0']

    def _gen_single_rnd_seq(self, seq_len, cutoff_length, scale, tries=10):
        """
        Internal function for random sequence generation.

        This function is intended for internal usage; for generating a single
        random sequence, use the standard `self._random_seq()` function
        with a `k` parameter of 1.
        """
        # Set the initial state for the random sequence, which is just the
        # pad symbol times the maximum preceding order.
        rnd_seq = (self._padsymbol,) * max(self._pre)

        # Repeatedly build a dictionary with the search space for next element
        # given the current value of `rnd_seq`, checking if we are able to
        # generate it (we might run into some unsolvable situation).
        gen_tries = 0
        while True:
            # Filter all the elements from `self._ngram_space` that match the
            # specified cutoff_length; as the conditional checking is a bit
            # expansive and not always used, we first check if we just make
            # a copy of the contents (casting the Counter to a dictionary,
            # which is fast) if no cutoff length is specified (i.e., if it
            # equals to 1).
            if cutoff_length == 1:
                sspace = dict(self._ngram_space)
            else:
                sspace = {key:value
                    for key, value in self._ngram_space.items()
                    if len(key) >= cutoff_length}

            # Filter all the items in the search space `sspace` that do not
            # start with what we have in the random sequence so far.
            sspace = {key:value for key, value in sspace.items()
                if key[:-1] == rnd_seq[-len(key)+1:]}

            # If the random sequence plus the new element would match the
            # sequence length, we can only use entries that end with the
            # padding symbol (the boundary symbol) if those are used
            # (max(self._post)>0); in all other cases, we must make sure that
            # the element is *not* a padding symbol. `seq_len` here already
            # includes the maximum padding on the left and at most one
            # padding on the right, if the model was trained as such.
            if max(self._post) > 0:
                if len(rnd_seq) + 1 == seq_len:
                    sspace = {key:value for key, value in sspace.items()
                        if key[-1] == self._padsymbol}
                else:
                    sspace = {key:value for key, value in sspace.items()
                        if key[-1] != self._padsymbol}

            # Scale the counts by key length (so we favor longer ngrams that
            # should be able to better capture the likelihood, while keeping
            # room for less frequent but observed elements) and by the
            # provided `scale`, which can be set to 1 for no effect.
            sspace = {key:(value*len(key))**scale
                for key, value in sspace.items()}

            # If we were unable to get a suitable searching space, the
            # generation failed and we must signal that; otherwise, let's
            # choose a random ngram for the search space and append its new
            # element (at the end, index -1) to our random sequence.
            if len(sspace) == 0:
                # We will get to this point if the generation failed because no
                # suitable search space was found. We just reset the random
                # sequence to the initial state and keep trying until we
                # exhaust the number of tries, retuning `None` on failure.
                # There would be better logics for this, but they might
                # impact the idea of having this function as a "primitive" if
                # we ever move to a different architecture.
                gen_tries += 1
                if gen_tries < tries:
                    rnd_seq = (self._padsymbol,) * max(self._pre)
                else:
                    return None
            else:
                pop = list(sspace.keys())
                w = list(sspace.values())
                rnd_seq += (_random_choices(pop, w)[0][-1], )

                # If we are now at the requested length for the random
                # sequence, let's exit the loop and the return the random
                # sequence.
                if len(rnd_seq) == seq_len:
                    return rnd_seq

    # TODO: what if unable to find words?
    # TODO: add seed
    def random_seq(self, scale=2, only_longest=False, k=1, tries_scale=10, seed=None):
        # Set the random seed.
        random.seed(seed)

        # Setup the cutoff length according to whether we should only use the
        # longest possible ngrams or not. The filtering is done inside the
        # main loop by selecting only ngrams which are equal in length or
        # larger than the specified `cutoff_length`. The unitary element
        # accounts for having at least the element being generated.
        if only_longest:
            cutoff_length = max(self._pre) + 1 + max(self._post)
        else:
            cutoff_length = 1

        # Initialize the list for holding the random sequences and cache
        # some values used repeatedly.
        rnd_seqs = []
        len_pop = list(self._seqlens.keys())
        len_w = list(self._seqlens.values())
        for i in range(k*tries_scale):
            # The length of the random sequence is selected in advance from
            # self._seqlens (plus the left and right orders, if any, for
            # dealing with padding symbols). As `_random_choices()` always
            # returns a list, we extract the first element.
            # NOTE: While it would be possible to populate a list of random
            #       lengths from _random_choices(), we use this less elegant
            #       solution of selecting a new random one each time to make
            #       sure that we are not stuck with a given length in case it
            #       is problematic given the ngram search space.
            rnd_seq_len = _random_choices(len_pop, len_w)[0]
            rnd_seq_len += max(self._pre)
            if max(self._post) > 0:
                rnd_seq_len += 1

            # Try to generate a random sequence and append it if successful.
            _rnd_seq = self._gen_single_rnd_seq(rnd_seq_len,
                                                cutoff_length, scale)
            if _rnd_seq:
                rnd_seqs.append(_rnd_seq)

            # Break the loop if we already got what we wanted.
            if len(rnd_seqs) == k:
                break

        # Return the randomly generated sequences, if any, without the
        # padding symbols; we don't need to query each element in each sequence
        # for identity and can just cut with the right indexes.
        return [rnd_seq[max(self._pre):-1] for rnd_seq in rnd_seqs]


# This method with zip, besides returning an iterator as desired, is faster
# than both the previous lingpy implementation and the one in NLTK; as this is
# the core of the ngram methods, it is important to have at least this
# primitive as fast as possible. This is intentionally not defaulting to any
# value for the order, so that users won't confuse a given order to all
# orders up to and including the given one.
def get_n_ngrams(sequence, order, pad_symbol=PAD_SYMBOL):
    """
    Build an iterator for collecting all ngrams of a given order.

    The sequence can optionally be padded with boundary symbols which are
    equal for before and and after sequence boundaries.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    order: int
        The order of the ngrams to be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    Returns
    -------
    out: iterable
        An iterable over the ngrams of the sequence, returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in get_n_ngrams(sent, 2):
    ...     print(ngram)
    ...
    ('$$$', 'Insurgents')
    ('Insurgents', 'killed')
    ('killed', 'in')
    ('in', 'ongoing')
    ('ongoing', 'fighting')
    ('fighting', '$$$')

    >>> for ngram in get_n_ngrams(sent, 1):
    ...     print(ngram)
    ...
    ('Insurgents',)
    ('killed',)
    ('in',)
    ('ongoing',)
    ('fighting',)

    >>> for ngram in get_n_ngrams(sent, 0):
    ...     print(ngram)
    ...
    """

    # Convert to a tuple, for faster computation, and pad the sequence if
    # needed. The test for `pad_symbol` is a bit more expansive (None and not
    # all False values) as we should allow False values to be padded for some
    # situations.
    seq = _seq_as_tuple(sequence)
    if pad_symbol is not None:
        seq = chain((pad_symbol,)* (order-1), seq, (pad_symbol,) * (order-1))
        seq = tuple(seq)

    # We generate the collection of ngrams for counting occurences with Python
    # `zip()` function, so we can rely on internal C-code for speeding things
    # up. What we do is build a list of arguments for the function by using a
    # list comprehension over variable `order` and, then, decompose such list
    # when passing it to the function. The list comprehension is built so that,
    # for a sequence such as the characters in "Markov" (here without
    # boundaries, but you should get the point) and an order of 3, we'll have:
    #   [['M', 'a', 'r', 'k', 'o', 'v'],
    #    ['a', 'r', 'k', 'o', 'v'],
    #    ['r', 'k', 'o', 'v'],
    #    ['k', 'o', 'v']]
    #
    # From which we zip all possible combinations.

    for ngram in zip(*[seq[i:] for i in range(order)]):
        yield ngram


def get_all_ngrams(sequence, orders=None, pad_symbol=PAD_SYMBOL):
    """
    Build an iterator for collecting all ngrams of a given set of orders.

    If no set of orders (i.e., "lengths") is provided, this will collect all
    possible ngrams in the sequence.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    orders: list
        An optional list of the orders of the ngrams to be collected. Can be
        larger than the length of the sequence, in which case the latter will
        be padded accordingly if requested. Defaults to the collection of all
        possible ngrams in the sequence with the minimum padding.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    Returns
    -------
    out: iterable
        An iterable over the ngrams of the sequence, returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents were killed"
    >>> for ngram in get_all_ngrams(sent):
    ...     print(ngram)
    ...
    ('Insurgents',)
    ('were',)
    ('killed',)
    ('$$$', 'Insurgents')
    ('Insurgents', 'were')
    ('were', 'killed')
    ('killed', '$$$')
    ('$$$', '$$$', 'Insurgents')
    ('$$$', 'Insurgents', 'were')
    ('Insurgents', 'were', 'killed')
    ('were', 'killed', '$$$')
    ('killed', '$$$', '$$$')
    """

    # Convert to a tuple, for faster computation, compute the orders (if they
    # were not given), and pad the sequence if requested.
    seq = _seq_as_tuple(sequence)
    if not orders:
        orders = range(len(seq)+1)

    for order in orders:
        for ngram in get_n_ngrams(seq, order, pad_symbol):
            yield ngram


def get_skipngrams(sequence, order, max_gaps, pad_symbol=PAD_SYMBOL, single_gap=True):
    """
    Build an iterator for collecting all skip ngrams of a given length.

    The function requires an information of the length of the skip ngrams to be
    collected, allowing to either collect ngrams with an unlimited number
    of gap openings (as described and implemented in Guthrie et al. 2006) or
    with at most one gap opening.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected. Must not include
        "None" as an element, as it is used as a sentinel during skip ngram
        collection following the implementation offered by Bird et al. 2018
        (NLTK), which is a de facto standard.

    order: int
        The order of the ngrams to be collected (parameter "n" in Guthrie et
        al. 2006).

    max_gaps: int
        The maximum number of gaps in the ngrams to be collected (parameter "k"
        in Guthrie et al. 2006).

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    single_gap: boolean
        An optional logic value indicating if multiple gap openings are to be
        allowed, as in Guthrie et al. (2006) and Bird et al. (2018), or if at
        most one gap_opening is to be allowed. Defaults to True.

    Returns
    -------
    out: iterable
        An iterable over the ngrams of the sequence, returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in get_skipngrams(sent, 2, 2):
    ...     print(ngram)
    ...
    ('$$$', 'Insurgents')
    ('Insurgents', 'killed')
    ('killed', 'in')
    ('in', 'ongoing')
    ('ongoing', 'fighting')
    ('fighting', '$$$')
    ('$$$', 'killed')
    ('Insurgents', 'in')
    ('killed', 'ongoing')
    ('in', 'fighting')
    ('ongoing', '$$$')
    ('$$$', 'in')
    ('Insurgents', 'ongoing')
    ('killed', 'fighting')
    ('in', '$$$')
    >>> for ngram in get_skipngrams(sent, 2, 2, single_gap=False):
    ...     print(ngram)
    ...
    ('$$$', 'Insurgents')
    ('$$$', 'killed')
    ('$$$', 'in')
    ('Insurgents', 'killed')
    ('Insurgents', 'in')
    ('Insurgents', 'ongoing')
    ('killed', 'in')
    ('killed', 'ongoing')
    ('killed', 'fighting')
    ('in', 'ongoing')
    ('in', 'fighting')
    ('in', '$$$')
    ('ongoing', 'fighting')
    ('ongoing', '$$$')
    ('fighting', '$$$')
    """

    # Check skip ngram length, which by definition must be at least two (one
    # element for the left side and one for the right one)
    if order < 2:
        raise ValueError("Skip ngram order must be at least 2.")

    # Convert to a tuple, if needed, and pad the sequence if requested while
    # caching the sequence length; please note that in this case we are
    # caching the sequence length *after* padding, so skip ngrams where one
    # of the sides is composed entirely of padded symbols will be included in
    # the collection. We don't do this with the primitive "get_n_ngrams()"
    # because we will later add our own temporary padding for ngram filtering;
    # this also ends up speeding things a little bit, as the conversion to
    # a tuple is only perfomed once.
    seq = _seq_as_tuple(sequence)
    if pad_symbol:
        seq = tuple(chain((pad_symbol,)* (order-1),
                          seq,
                          (pad_symbol,) * (order-1)))
    len_seq = len(seq)

    # The logic for obtaining skip ngrams is different if we allow for multiple
    # gaps (as proposed and implemented by both Guthrie et al. 2006 and Bird et
    # al. in NLTK, whose code logic in module nltk.util at
    # http://www.nltk.org/_modules/nltk/util.html is closely followed here), or
    # if we allow for a single gap, especially considering that we cannot
    # collet repeated ngrams (with a gap of zero, an ngram for preceding length
    # 1 and following length 2 is equal to an ngram of preceding length 2 and
    # following length 1).
    if not single_gap:
        # We pad the `sequence` with None symbols to the right, so we can
        # filter during the list comprehension. Please note that this is *not*
        # the user-requested padding, but an internal and inexpansive way to
        # account for the end of the ngram; also note that this will fail if
        # the sequence itself holds None symbols.
        # NOTE: To solve the problems when/if the sequence itself holds None
        #       symbols, we could translate Nones to a tempory value and remap
        #       it when yielding; while expansive, this is still more effective
        #       (and clear) than other solutions which would not allow using
        #       `all()` for most sequence computations.
        _temp = chain(seq, (None,) * order)
        for ngram in get_n_ngrams(_temp, order+max_gaps, pad_symbol=None):
            head = ngram[:1] # cache for all combinations
            for comb in [tail for tail in combinations(ngram[1:], order-1) if all(tail)]:
                yield head + comb
    else:
        # Iterate over all the possible gap lengths, including length zero.
        # Length zero requires a different logic (it is actually just returning
        # all subsequent ngrams) in order not to yield the same subsequence
        # (for example, for length 1+2 and 2+1, which are obviously identical).
        # One alternative would be to add a function paramenter allowing the
        # user to circumvent this behaviour, but no situation where this would
        # be required or at least suggested seem to exist (we would be
        # distributing more probability to ngrams with no gaps, which is
        # exactly what skip ngrams are intended to correct).
        for gap_width in range(max_gaps+1):
            if gap_width == 0:
                # No pad-symbol here, as the sequence was padded before.
                for ngram in get_n_ngrams(seq, order, pad_symbol=None):
                    yield ngram
            else:
                # We iterate over all possible left and right lengths,
                # making sure we always have at least one on each side.
                for left_width in range(1, order):
                    # Build the pattern we are extracting
                    pattern = (True,) * left_width + \
                              (False,) * gap_width + \
                              (True,) * (order - left_width) # right width

                    # Iterate over the sequence while applying the pattern;
                    # `idx` is the starting index in `sequence`, `skip` the
                    # element index in `pattern`
                    for idx in range(len_seq-order-gap_width+1):
                        yield tuple(seq[idx+skip] for skip, keep in enumerate(pattern) if keep)

def get_posngrams(sequence, pre_order=0, post_order=0, pad_symbol=PAD_SYMBOL, elm_symbol=ELM_SYMBOL):
    """
    Build an iterator for collecting all positional ngrams of a sequence.

    The preceding and a following orders (i.e., "contexts") must always be
    informed. The elements of the iterator include a tuple of the context,
    which can be hashed as any tuple, the transition symbol, and the position
    of the symbol in the sequence. Such output is primarily intended for
    state-by-state relative likelihood computations with stochastics models.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    pre-order: int
        An optional integer specifying the length of the preceding context.
        Defaults to zero.

    post-order: int
        An optional integer specifying the length of the following context.
        Defaults to zero.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    elm_symbol: object
        An optional symbol to be used as transition symbol replacement in the
        context tuples (the first element in the returned iterator). Defaults
        to "###".

    Returns
    -------
    out: iterable
        An iterable over the positional ngrams of the sequence, returned as
        tuples whose elements are: (1) a tuple representing the context
        (thus including preceding context, the transition symbol, and the
        following context), (2) an object with the value of the transition
        symbol, and (3) the index of the transition symbol in the sequence.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in get_posngrams(sent, 2, 1):
    ...     print(ngram)
    ...
    (('$$$', '$$$', '###', 'killed'), 'Insurgents', 0)
    (('$$$', 'Insurgents', '###', 'in'), 'killed', 1)
    (('Insurgents', 'killed', '###', 'ongoing'), 'in', 2)
    (('killed', 'in', '###', 'fighting'), 'ongoing', 3)
    (('in', 'ongoing', '###', '$$$'), 'fighting', 4)
    """

    # Cache the complexive order for the ngram from the sum of the pre- and
    # post- orders (with an additional one, the state under actual
    # observation).
    order = pre_order + 1 + post_order

    # Pad the sequence if requested and cache the sequence length for later
    # deciding whether to include an ngram based on the state index (thus
    # excluding ngrams centered in padded symbols, which would otherwise be
    # impossible to identify). Please note that in this case of positional
    # ngrams (unlike skip ngrams, for example), the sequence length is cached
    # *before* padding precisely in order to allow the filtering of elements.
    seq = _seq_as_tuple(sequence)
    if pad_symbol:
        seq = chain((pad_symbol,)* pre_order, seq, (pad_symbol,) * post_order)
        seq = tuple(seq)

    # We obtain all the subsequences of the order we desire by asking for all
    # the ngrams of the given order when the sequence is not addionally padded
    # (of course, it will already have been padded, if the user so requested,
    # by this time).
    subseqs = get_n_ngrams(seq, order, pad_symbol=None)

    # We can now collect all the skipping sequences, caching the various
    # indexes for quicker extraction. We chain from iterables as this is
    # faster for such a primitive function.
    elem_idx = -1 - post_order
    postctx_idx = order - post_order
    for state_idx, subseq in enumerate(subseqs):
        yield (
            # pre-context + element + post_context
            (subseq[:elem_idx] + (elm_symbol,) + subseq[postctx_idx:]),
            # transition
            subseq[elem_idx],
            # state index
            state_idx)

def get_all_posngrams(sequence, pre_orders, post_orders, pad_symbol=PAD_SYMBOL, elm_symbol=ELM_SYMBOL):
    """
    Build an iterator for collecting all positional ngrams of a sequence.

    The elements of the iterator, as returned by "get_posngrams()", include a
    tuple of the context, which can be hashed (as any tuple), the transition
    symbol, and the position of the symbol in the sequence. Such output is
    primarily intended for state-by-state relative likelihood computations with
    stochastics models, and can be approximated to a collection of "shingles".

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    pre-orders: int or list
        An integer with the maximum length of the preceding context or a list
        with all preceding context lengths to be collected. If an integer is
        passed, all lengths from zero to the informed one will be collected.

    post-orders: int or list
        An integer with the maximum length of the following context or a list
        with all following context lengths to be collected. If an integer is
        passed, all lengths from zero to the informed one will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    elm_symbol: object
        An optional symbol to be used as transition symbol replacement in the
        context tuples (the first element in the returned iterator). Defaults
        to "###".

    Returns
    -------
    out: iterable
        An iterable over the positional ngrams of the sequence, returned as
        tuples whose elements are: (1) a tuple representing the context (thus
        including preceding context, the transition symbol, and the following
        context), (2) an object with the value of the transition symbol, and
        (3) the index of the transition symbol in the sequence.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents were killed"
    >>> for ngram in get_all_posngrams(sent, 2, 1):
    ...     print(ngram)
    ...
    (('###',), 'Insurgents', 0)
    (('###',), 'were', 1)
    (('###',), 'killed', 2)
    (('###', 'were'), 'Insurgents', 0)
    (('###', 'killed'), 'were', 1)
    (('###', '$$$'), 'killed', 2)
    (('$$$', '###'), 'Insurgents', 0)
    (('Insurgents', '###'), 'were', 1)
    (('were', '###'), 'killed', 2)
    (('$$$', '###', 'were'), 'Insurgents', 0)
    (('Insurgents', '###', 'killed'), 'were', 1)
    (('were', '###', '$$$'), 'killed', 2)
    (('$$$', '$$$', '###'), 'Insurgents', 0)
    (('$$$', 'Insurgents', '###'), 'were', 1)
    (('Insurgents', 'were', '###'), 'killed', 2)
    (('$$$', '$$$', '###', 'were'), 'Insurgents', 0)
    (('$$$', 'Insurgents', '###', 'killed'), 'were', 1)
    (('Insurgents', 'were', '###', '$$$'), 'killed', 2)
    """

    # We don't need to convert `sequence` into a tuple or pad it here, as this
    # will be performed by `get_posngrams()`. While we could do this in advance
    # and cache the results, this complicates things a bit and a quick
    # experimentation showed no real improvement in perfomance, even when
    # simulating with large datasets (we'd still need to perform a conditional
    # check on the sequence type in order to profit from the cache, which is
    # expensive and is otherwise performed internally by C-code).

    # For both pre- and post-context, we will interact over all lengths if we
    # receive a list, or build a range of such lengths if an integer is
    # received (for this reason, we add a unit to the range, so that the top
    # value when passing an integer is compatible to max() when passing a list).
    if isinstance(pre_orders, int):
        pre_orders = range(pre_orders + 1)
    if isinstance(post_orders, int):
        post_orders = range(post_orders + 1)

    # Collect all ngrams...
    ngrams = [
        get_posngrams(sequence, pre_order, post_order, pad_symbol, elm_symbol)
        for pre_order, post_order
        in product(pre_orders, post_orders)]

    # ...and yield them; there is probably a way of having this a bit
    # more functional even in Python, but it would likely complicate the code
    # too much and unnecessarily.
    for ngram in chain.from_iterable(ngrams):
        yield ngram


# Define partial functions

bigrams = partial(get_n_ngrams, order=2)
bigrams.__doc__ = """
    Build an iterator for collecting all bigrams of a sequence.

    The sequence is padded by default.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the bigrams will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    Returns
    -------
    out: iterable
        An iterable over the bigrams of the sequence, returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in bigrams(sent):
    ...     print(ngram)
    ...
    ('$$$', 'Insurgents')
    ('Insurgents', 'killed')
    ('killed', 'in')
    ('in', 'ongoing')
    ('ongoing', 'fighting')
    ('fighting', '$$$')
    """


trigrams = partial(get_n_ngrams, order=3)
trigrams.__doc__ = """
    Build an iterator for collecting all trigrams of a  sequence.

    The sequence is padded by default.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the trigrams will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    Returns
    -------
    out: iterable
        An iterable over the trigrams of the sequence, returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in trigrams(sent):
    ...     print(ngram)
    ...
    ('$$$', '$$$', 'Insurgents')
    ('$$$', 'Insurgents', 'killed')
    ('Insurgents', 'killed', 'in')
    ('killed', 'in', 'ongoing')
    ('in', 'ongoing', 'fighting')
    ('ongoing', 'fighting', '$$$')
    ('fighting', '$$$', '$$$')
    """


fourgrams = partial(get_n_ngrams, order=4)
fourgrams.__doc__ = """
    Build an iterator for collecting all fourgrams of a sequence.

    The sequence is padded by default.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the fourgrams will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".

    Returns
    -------
    out: iterable
        An iterable over the fourgrams of the sequence, returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in fourgrams(sent):
    ...     print(ngram)
    ...
    ('$$$', '$$$', '$$$', 'Insurgents')
    ('$$$', '$$$', 'Insurgents', 'killed')
    ('$$$', 'Insurgents', 'killed', 'in')
    ('Insurgents', 'killed', 'in', 'ongoing')
    ('killed', 'in', 'ongoing', 'fighting')
    ('in', 'ongoing', 'fighting', '$$$')
    ('ongoing', 'fighting', '$$$', '$$$')
    ('fighting', '$$$', '$$$', '$$$')
    """
