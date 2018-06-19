# *-* coding: utf-8 *-*
"""
This module provides various methods for generating and
collecting n-grams from sequences.

"""
# TODO: write above ngrams, pngrams, skip ngrams


# TODO: check the import in __init__.py

from itertools import chain, combinations, product

_ELEMENT = "###"

# TODO: remove from sound_classes.py later
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
# important to have at least this as fast as possible.
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

    return zip(*[seq[i:] for i in range(order)])
    

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

    """
    
    if not orders:
        orders = range(len(sequence)+1)
        
    for order in orders:
        for ngram in get_n_ngrams(sequence, order, pad_symbol):
            yield ngram


def get_pgrams(sequence, pre_order=0, post_order=0, pad_symbol=None):
    """
    Collect a single pre- and post- order ngram set from a sequence.
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
    len_seq = len(seq)
    if pad_symbol:
        seq = chain((pad_symbol,)* pre_order, seq, (pad_symbol,) * post_order)
        seq = tuple(seq)
        
    # We obtain all the subsequences of the order we desire by
    # asking for the all the ngrams of the given order when the
    # sequence is not addionally padded (of course, it will already
    # have been padded, if the user so requested, by this time).
    subseqs = get_n_ngrams(seq, order, pad_symbol=None)

#    print(len(list(subseqs)))
#    print(list(subseqs))

    # We can now collect all the skipping sequences, caching the
    # various indexes for quicker extraction. Given the way
    # Python indexes lists, we need to perform a conditional
    # extraction based on the value of post_order (otherwise
    # we would get a `postctx_idx` of -1+1, i.e. zero, which
    # would return the entire subsequence; work-arounds would
    # make the code too complex for something so simple, and
    # this has the advantage of being faster in cases of
    # only subsequent ngrams, which are more common).
    if not post_order:
        ngrams = (
            (tuple(subseq[:-1] + (_ELEMENT,)), subseq[-1], state_idx + pre_order)
            for state_idx, subseq
            in enumerate(subseqs))
    else:
        elem_idx = -1 - post_order
        prectx_idx = elem_idx - pre_order
        postctx_idx = -post_order
        ngrams = (
            # pre-context + element + post_context
            (tuple(subseq[prectx_idx:elem_idx] + (_ELEMENT,) + subseq[postctx_idx:]),
             # transition
             subseq[elem_idx],
             # state index, discounting the start boundaries
             state_idx + pre_order)
            # loop
            for state_idx, subseq in enumerate(subseqs)
        )

    return ngrams


def collect_ngrams_multiorder(seq, pre_order, post_order, pad_symbol="$$$"):
    """
    Collect various order ngrams from a sequence.
    """
    # For both pre- and post-context, we will interact over all
    # lengths if we receive a list, or build a range of such lengths
    # if an integer is received (for this reason, we add a unit to
    # the range, so that the top value when passing an integer is
    # compatible to max() when passing a list).
    if isinstance(pre_order, int):
        pre_order = range(pre_order + 1)
    if isinstance(post_order, int):
        post_order = range(post_order + 1)

    # Collect all ngrams...
    ngrams = [collect_ngrams(seq, pre_length, post_length, pad_symbol)
              for pre_length, post_length
              in product(pre_order, post_order)]

    # ...and flatten the list, so it can easily be passed to a Counter later
    # (and so that it is compatible with what is returned by
    # `collect_ngrams()`)
    ng = [ngram for order_ngrams in ngrams for ngram in order_ngrams]

    return ng


def skip_ngrams(sequence, n, k, pad_symbol=None, subsequent=True):
    # Check skip ngram length, which by definition must be at
    # least two (one element for the left side and one for the
    # right one)
    if n < 2:
        raise ValueError("Skip ngram order must be at least 2.")

    # Pad the sequence if requested and cache the sequence length;
    # please note that in this case we are caching the
    # sequence length *after* padding, so skip ngrams where one
    # of the sides is composed entirely of padded symbols
    # will be included.
    if pad_symbol:
        sequence = chain((pad_symbol,)* (n-1), sequence, (pad_symbol,) * (n-1))
        sequence = list(sequence)
    len_seq = len(sequence)
        
    # The logic for obtaining skip ngrams is different if we allow
    # for multiple gaps (as proposed by Guthrie et al. 2006 and
    # implemented, among others, by Bird et al. 2004, whose
    # implementation is followed here), or if we allow
    # for a single gap, especially considering that we should not
    # offer repeated ngrams. 
    if not subsequent:
        # We pad the `sequence` with None symbols to the right,
        # so we can filter during the list comprehension.
        # Please note that this is *not* the user-requested
        # padding, but an internal and inexpansive way to
        # account for the end of the ngram;
        # please note that this will fail if the sequence itself
        # holds None symbols.
        # TODO: To solve the problems when/if the sequence itself
        #       holds None symbols, we could translate Nones to
        #       a tempory value and remap it when yielding; while
        #       expansive, this is still more effective than
        #       other solutions which would not allow using
        #       `all()` for most sequence computations.
        _temp = chain(sequence, (None,) * n)
        for ngram in get_ngrams(_temp, n+k, pad_symbol=None):
            head = ngram[:1]
            combs = [tail for tail in combinations(ngram[1:], n-1) if all(tail)]
            for comb in combs:
                yield head + comb
    else:
        # Iterate over all the possible gap lengths, including
        # length zero. Length zero requires a different logic (it is
        # actually just returning all subsequent ngrams) in order
        # not to yield the same subsequence (for example, for length
        # 1+2 and 2+1, which are obviously identical). One alternative,
        # if needed, is to add a function parameter circumventing this.
        for gap_width in range(0, k+1):
            if gap_width == 0:
                for ngram in get_ngrams(sequence, n, pad_symbol=None):
                    yield ngram
            else:
                # We iterate over all possible left and right lengths,
                # making sure we always have at least one on each side.
                for left_width in range(1, n):
                    # Compute the right width
                    right_width = n - left_width
                    # Build the pattern we are extracting
                    pattern = (True,) * left_width + \
                              (False,) * gap_width + \
                              (True,) * right_width
                    # Iterate over the sequence while applying the pattern;
                    # `i` is the starting index in `sequence`, `j` the
                    # element index in `pattern`
                    for i in range(len_seq-len(pattern)+1):
                        ngram = [sequence[i+j] for j, keep in enumerate(pattern) if keep]
                        yield tuple(ngram)
