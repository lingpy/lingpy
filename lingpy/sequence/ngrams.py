# *-* coding: utf-8 *-*
"""
This module provides various methods for generating and
collecting n-grams from sequences. Standard ngrams,
skip ngrams, and positional ngrams can be collected.
"""

from itertools import chain, combinations, product

def _seq_as_tuple(sequence):
    """
    Internal function for automatically converting a
    string sequence to a tuple, if needed.

    Parameters
    ----------
    sequence: list or str
        The sequence that shall be converted into an
        iterable.

    Returns
    -------
    out: tuple
        A tuple of the sequence.
    """

    # We first check for datatype and then for a space,
    # as the first test is faster (and evaluation is lazy).
    if isinstance(sequence, str) and ' ' in sequence:
        return tuple(sequence.split(' '))

    return tuple(sequence)

# This method with zip, besides returning an iterator as desired,
# is faster than both the previous lingpy implementation and the
# one in NLTK; as this is the core of the ngram methods, it is
# important to have at least this primitive as fast as possible.
# This is intentionally not defaulting to any value for the
# order, so that users won't confuse a given order to all
# orders up to and including the given one.
def get_n_ngrams(sequence, order, pad_symbol='$'):
    """
    Build an iterator for collecting all ngrams of a given
    order from a sequence. The sequence can optionally be padded
    with boundary symbols which are equal for before and and
    after sequence boundaries.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    order: int
        The order of the ngrams to be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    Returns
    -------
    out: iterable
        An iterable over the ngrams of the sequence,
        returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in get_n_ngrams(sent, 2):
    ...     print(ngram)
    ...
    ('$', 'Insurgents')
    ('Insurgents', 'killed')
    ('killed', 'in')
    ('in', 'ongoing')
    ('ongoing', 'fighting')
    ('fighting', '$')

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

    # Convert to a tuple, for faster computation, and pad the
    # sequence if needed. The test for `pad_symbol` is a bit
    # more expansive (None and not all False values) as we
    # should allow False values to be padded for some situations.
    seq = _seq_as_tuple(sequence)
    if pad_symbol is not None:
        seq = chain((pad_symbol,)* (order-1), seq, (pad_symbol,) * (order-1))
        seq = tuple(seq)

    # We generate the collection of ngrams for counting occurences
    # with Python `zip()` function, so we can
    # rely on internal C-code for speeding things up. What we do
    # is build a list of arguments for the function by using a
    # list comprehension over variable `order` and, then,
    # decompose such list when passing it to the function.
    # The list comprehension is built so that, for
    # a sequence such as the characters in "Markov" (here without
    # boundaries, but you should get the point) and an order of 3,
    # we'll have:
    #   [['M', 'a', 'r', 'k', 'o', 'v'],
    #    ['a', 'r', 'k', 'o', 'v'],
    #    ['r', 'k', 'o', 'v'],
    #    ['k', 'o', 'v']]
    #
    # From which we zip all possible combinations.

    for ngram in zip(*[seq[i:] for i in range(order)]):
        yield ngram

def bigrams(sequence, pad_symbol='$'):
    """
    Build an iterator for collecting all bigrams of a
    sequence, padding the sequence by default.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the bigrams will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    Returns
    -------
    out: iterable
        An iterable over the bigrams of the sequence,
        returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in bigrams(sent):
    ...     print(ngram)
    ...
    ('$', 'Insurgents')
    ('Insurgents', 'killed')
    ('killed', 'in')
    ('in', 'ongoing')
    ('ongoing', 'fighting')
    ('fighting', '$')
    """

    for bigram in get_n_ngrams(sequence, 2, pad_symbol):
        yield bigram

def trigrams(sequence, pad_symbol='$'):
    """
    Build an iterator for collecting all trigrams of a
    sequence, padding the sequence by default.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the trigrams will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    Returns
    -------
    out: iterable
        An iterable over the trigrams of the sequence,
        returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in trigrams(sent):
    ...     print(ngram)
    ...
    ('$', '$', 'Insurgents')
    ('$', 'Insurgents', 'killed')
    ('Insurgents', 'killed', 'in')
    ('killed', 'in', 'ongoing')
    ('in', 'ongoing', 'fighting')
    ('ongoing', 'fighting', '$')
    ('fighting', '$', '$')
    """

    for trigram in get_n_ngrams(sequence, 3, pad_symbol):
        yield trigram


def fourgrams(sequence, pad_symbol='$'):
    """
    Build an iterator for collecting all fourgrams of a
    sequence, padding the sequence by default.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the fourgrams will be collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    Returns
    -------
    out: iterable
        An iterable over the fourgrams of the sequence,
        returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in fourgrams(sent):
    ...     print(ngram)
    ...
    ('$', '$', '$', 'Insurgents')
    ('$', '$', 'Insurgents', 'killed')
    ('$', 'Insurgents', 'killed', 'in')
    ('Insurgents', 'killed', 'in', 'ongoing')
    ('killed', 'in', 'ongoing', 'fighting')
    ('in', 'ongoing', 'fighting', '$')
    ('ongoing', 'fighting', '$', '$')
    ('fighting', '$', '$', '$')
    """

    for fourgram in get_n_ngrams(sequence, 4, pad_symbol):
        yield fourgram


def get_all_ngrams(sequence, orders=None, pad_symbol='$'):
    """
    Build an iterator for collecting all ngrams of a sequence
    of a given set of orders (i.e., "lengths"). If no set of
    orders is provided, this will collect all possible ngrams
    in the sequence.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    orders: list
        An optional list of the orders of the ngrams to
        be collected. Can be larger than the length of the
        sequence, in which case the latter will be padded
        accordingly if requested. Defaults to the collection
        of all possible ngrams in the sequence with the
        minimum padding.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    Returns
    -------
    out: iterable
        An iterable over the ngrams of the sequence,
        returned as tuples.

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
    ('$', 'Insurgents')
    ('Insurgents', 'were')
    ('were', 'killed')
    ('killed', '$')
    ('$', '$', 'Insurgents')
    ('$', 'Insurgents', 'were')
    ('Insurgents', 'were', 'killed')
    ('were', 'killed', '$')
    ('killed', '$', '$')
    """

    # Convert to a tuple, for faster computation, compute the
    # orders (if they were not given), and pad the sequence if
    # requested.
    seq = _seq_as_tuple(sequence)
    if not orders:
        orders = range(len(seq)+1)

    for order in orders:
        for ngram in get_n_ngrams(seq, order, pad_symbol):
            yield ngram


def get_skipngrams(sequence, order, max_gaps, pad_symbol='$', single_gap_opening=True):
    """
    Build an iterator for collecting all skip ngrams of a sequence
    of a given length and of a maximum number of gaps,
    with with unlimited number of gap openings (as described in
    Guthrie et al. 2006) or with at most one gap opening.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.
        Must not include "None" as an element, as it is used
        a sentinel during skip ngram collection following the
        implementation in Bird et al. 2018 (NLTK) which is
        a de facto standard.

    order: int
        The order of the ngrams to be collected (parameter
        "n" in Guthrie et al. 2006).

    max_gaps: int
        The maximum number of gaps in the ngrams to be
        collected (parameter "k" in Guthrie et al. 2006).

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    single_gap_opening: boolean
        An optional logic value indicating if multiple
        gap openings are to be allowed, as in Guthrie et
        al. (2006) and Bird et al. (2018), or if at
        most one gap_opening is to be allowed. Defaults
        to True.

    Returns
    -------
    out: iterable
        An iterable over the ngrams of the sequence,
        returned as tuples.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in get_skipngrams(sent, 2, 2):
    ...     print(ngram)
    ...
    ('$', 'Insurgents')
    ('Insurgents', 'killed')
    ('killed', 'in')
    ('in', 'ongoing')
    ('ongoing', 'fighting')
    ('fighting', '$')
    ('$', 'killed')
    ('Insurgents', 'in')
    ('killed', 'ongoing')
    ('in', 'fighting')
    ('ongoing', '$')
    ('$', 'in')
    ('Insurgents', 'ongoing')
    ('killed', 'fighting')
    ('in', '$')
    >>> for ngram in get_skipngrams(sent, 2, 2, single_gap_opening=False):
    ...     print(ngram)
    ...
    ('$', 'Insurgents')
    ('$', 'killed')
    ('$', 'in')
    ('Insurgents', 'killed')
    ('Insurgents', 'in')
    ('Insurgents', 'ongoing')
    ('killed', 'in')
    ('killed', 'ongoing')
    ('killed', 'fighting')
    ('in', 'ongoing')
    ('in', 'fighting')
    ('in', '$')
    ('ongoing', 'fighting')
    ('ongoing', '$')
    ('fighting', '$')
    """

    # Check skip ngram length, which by definition must be at
    # least two (one element for the left side and one for the
    # right one)
    if order < 2:
        raise ValueError("Skip ngram order must be at least 2.")

    # Convert to a tuple, if needed, and
    # pad the sequence if requested and cache the sequence length;
    # please note that in this case we are caching the
    # sequence length *after* padding, so skip ngrams where one
    # of the sides is composed entirely of padded symbols
    # will be included in the collection. We don't do this with
    # the primitive "get_n_ngrams()" because we will later add
    # our own temporary padding for ngram filtering; this also
    # ends up speeding things a little bit, as the conversion to
    # a tuple is only perfomed once.
    seq = _seq_as_tuple(sequence)
    if pad_symbol:
        seq = tuple(chain((pad_symbol,)* (order-1), seq, (pad_symbol,) * (order-1)))
    len_seq = len(seq)

    # The logic for obtaining skip ngrams is different if we allow
    # for multiple gaps (as proposed and implemented by both
    # Guthrie et al. 2006 and Bird et al. in NLTK, whose
    # code in module nltk.util at http://www.nltk.org/_modules/nltk/util.html
    # is closely followed here), or if we allow
    # for a single gap, especially considering that we cannot
    # collet repeated ngrams (with a gap of zero, an ngram for
    # preceding length 1 and following length 2 is equal to
    # an ngram of preceding length 2 and following length 1).
    if not single_gap_opening:
        # We pad the `sequence` with None symbols to the right,
        # so we can filter during the list comprehension.
        # Please note that this is *not* the user-requested
        # padding, but an internal and inexpansive way to
        # account for the end of the ngram; also
        # please note that this will fail if the sequence itself
        # holds None symbols.
        # NOTE: To solve the problems when/if the sequence itself
        #       holds None symbols, we could translate Nones to
        #       a tempory value and remap it when yielding; while
        #       expansive, this is still more effective than
        #       other solutions which would not allow using
        #       `all()` for most sequence computations.
        _temp = chain(seq, (None,) * order)
        for ngram in get_n_ngrams(_temp, order+max_gaps, pad_symbol=None):
            head = ngram[:1] # cache for all combinations
            for comb in [tail for tail in combinations(ngram[1:], order-1) if all(tail)]:
                yield head + comb
    else:
        # Iterate over all the possible gap lengths, including
        # length zero. Length zero requires a different logic (it is
        # actually just returning all subsequent ngrams) in order
        # not to yield the same subsequence (for example, for length
        # 1+2 and 2+1, which are obviously identical). One alternative,
        # would be to add a function paramenter allowing the user to
        # circumvent this behaviour, but no situation where this would
        # be required or at least suggested seem to exist (we would be
        # distributing more probability to ngrams with no gaps, which
        # is exactly what skip ngrams are intended to correct).
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

def get_posngrams(sequence, pre_order=0, post_order=0, pad_symbol='$', elm_symbol='###'):
    """
    Build an iterator for collecting all positional ngrams of a sequence
    of a given a preceding and a following orders (i.e., "contexts").
    The elements of the iterator include a tuple of the context, which
    can be hashed as any tuple, the transition symbol, and the
    position of the symbol in the sequence. Such output is
    primarily intended for state-by-state relative likelihood
    computations with stochastics models.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    pre-order: int
        An optional integer specifying the length of the
        preceding context. Default to zero.

    post-order: int
        An optional integer specifying the length of the
        following context. Default to zero.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    elm_symbol: object
        An optional symbol to be used as transition
        symbol replacement in the context tuples
        (the first element in the returned iterator).
        Defaults to "###".

    Returns
    -------
    out: iterable
        An iterable over the positional ngrams of the
        sequence, returned as tuples whose elements are:
        (1) a tuple with representing the context (thus
        including preceding context, the transition
        symbol, and the following context), (2) an
        object with the value of the transition symbol,
        and (3) the index of the transition symbol in
        the sequence.

    Examples
    --------
    >>> from lingpy.sequence import *
    >>> sent = "Insurgents killed in ongoing fighting"
    >>> for ngram in get_posngrams(sent, 2, 1):
    ...     print(ngram)
    ...
    (('$', '$', '###', 'killed'), 'Insurgents', 0)
    (('$', 'Insurgents', '###', 'in'), 'killed', 1)
    (('Insurgents', 'killed', '###', 'ongoing'), 'in', 2)
    (('killed', 'in', '###', 'fighting'), 'ongoing', 3)
    (('in', 'ongoing', '###', '$'), 'fighting', 4)
    """

    # Cache the complexive order for the ngram from the sum of the
    # pre- and post- orders (with an additional one, the state
    # under actual observation).
    order = pre_order + 1 + post_order

    # Pad the sequence if requested and cache the sequence length
    # for later deciding whether to include an ngram based on
    # the state index (thus excluding ngrams centered
    # in padded symbols, which would otherwise be impossible to
    # identify). Please note that in this case of positional
    # ngrams (unlike skip ngrams, for example), the sequence length
    # is cache *before* padding precisely in order to allow
    # the filtering of elements.
    seq = _seq_as_tuple(sequence)
    if pad_symbol:
        seq = chain((pad_symbol,)* pre_order, seq, (pad_symbol,) * post_order)
        seq = tuple(seq)

    # We obtain all the subsequences of the order we desire by
    # asking for the all the ngrams of the given order when the
    # sequence is not addionally padded (of course, it will already
    # have been padded, if the user so requested, by this time).
    subseqs = get_n_ngrams(seq, order, pad_symbol=None)

    # We can now collect all the skipping sequences, caching the
    # various indexes for quicker extraction. We chain from
    # iterables as this is faster for such a primitive function.
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

def get_all_posngrams(sequence, pre_orders, post_orders, pad_symbol='$', elm_symbol='###'):
    """
    Build an iterator for collecting all positional ngrams of a sequence
    of a given set of preceding and following orders (i.e., "contexts").
    The elements of the iterator, as returned by "get_posngrams()",
    include a tuple of the context, which
    can be hashed as any tuple, the transition symbol, and the
    position of the symbol in the sequence. Such output is
    primarily intended for state-by-state relative likelihood
    computations with stochastics models.

    Parameters
    ----------
    sequence: list or str
        The sequence from which the ngrams will be collected.

    pre-orders: int
        An integer with the maximum length of the preceding
        context or a list with all preceding context lengths
        to be collected. If an integer is passed, all
        lengths from zero to the informed one will be
        collected.

    post-orders: int
        An integer with the maximum length of the following
        context or a list with all preceding context lengths
        to be collected. If an integer is passed, all
        lengths from zero to the informed one will be
        collected.

    pad_symbol: object
        An optional symbol to be used as start-of- and
        end-of-sequence boundaries. The same symbol
        is used for both boundaries. Must be a value
        different from None, defaults to "$".

    elm_symbol: object
        An optional symbol to be used as transition
        symbol replacement in the context tuples
        (the first element in the returned iterator).
        Defaults to "###".

    Returns
    -------
    out: iterable
        An iterable over the positional ngrams of the
        sequence, returned as tuples whose elements are:
        (1) a tuple with representing the context (thus
        including preceding context, the transition
        symbol, and the following context), (2) an
        object with the value of the transition symbol,
        and (3) the index of the transition symbol in
        the sequence.

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
    (('###', '$'), 'killed', 2)
    (('$', '###'), 'Insurgents', 0)
    (('Insurgents', '###'), 'were', 1)
    (('were', '###'), 'killed', 2)
    (('$', '###', 'were'), 'Insurgents', 0)
    (('Insurgents', '###', 'killed'), 'were', 1)
    (('were', '###', '$'), 'killed', 2)
    (('$', '$', '###'), 'Insurgents', 0)
    (('$', 'Insurgents', '###'), 'were', 1)
    (('Insurgents', 'were', '###'), 'killed', 2)
    (('$', '$', '###', 'were'), 'Insurgents', 0)
    (('$', 'Insurgents', '###', 'killed'), 'were', 1)
    (('Insurgents', 'were', '###', '$'), 'killed', 2)
    """

    # We don't need to convert `sequence` into a tuple or pad it
    # here, as this will be performed by `get_posngrams()`.
    # While we could do this in advance and cache the results,
    # this complicates things a bit and a quick experimentation
    # showed no real improvement in perfomance, even when
    # simulating with large datasets (we'd still need to
    # perform a conditional check on the sequence type in order
    # to profit from the cache, which is expensive and is
    # otherwise performed internally by C-code).

    # For both pre- and post-context, we will interact over all
    # lengths if we receive a list, or build a range of such lengths
    # if an integer is received (for this reason, we add a unit to
    # the range, so that the top value when passing an integer is
    # compatible to max() when passing a list).
    if isinstance(pre_orders, int):
        pre_orders = range(pre_orders + 1)
    if isinstance(post_orders, int):
        post_orders = range(post_orders + 1)

    # Collect all ngrams...
    ngrams = [get_posngrams(sequence, pre_order, post_order, pad_symbol, elm_symbol)
              for pre_order, post_order
              in product(pre_orders, post_orders)]

    # ...and yield them; there is probably a way of having this a bit
    # more functional even in Python, but it would likely complicate the code
    # too much and unnecessarily.
    for ngram in chain.from_iterable(ngrams):
        yield ngram
