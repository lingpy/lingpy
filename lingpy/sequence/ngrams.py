from itertools import combinations, chain, zip_longest, product

def get_ngrams(sequence, n, pad_symbol='$'):
    sequence = iter(sequence)

    if pad_symbol:
        sequence = chain((pad_symbol,)* (n-1), sequence, (pad_symbol,) * (n-1))

    history = [next(sequence) for i in range(n-1)]

    for item in sequence:
        history.append(item)
        yield tuple(history)
        del history[0]

def all_grams(sequence, min_len=1, max_len=None, pad_symbol='$'):
    if not max_len:
        max_len = len(sequence)

    for n in range(min_len, max_len+1):
        for ngram in get_ngrams(sequence, n, pad_symbol):   
            yield ngram


####################################################################

# We define some functions for collecting substrings, here called n-grams even though they are not
# necessarily proper ngrams. In all cases, we specify a pre- and a post- order for the amount of
# information before and after each state to be collected.

_ELEMENT = "###"
def collect_ngrams(seq, pre_order=0, post_order=0, pad_symbol=None):
    """
    Collect a single pre- and post- order ngram set from a sequence.
    """
    
    # Cache the sequence length, so we can later decided whether to
    # include an ngram based on the state index, and cache complexive
    # order for the ngram from the sum of the pre- and post- orders
    # (with an additional element, the referential one).
    len_seq = len(seq)
    order = pre_order + post_order + 1

    # Pad the sequence if requested; please note that in this case
    # (unlike skip ngram computation, for example), the sequence
    # length is cached *before* the paddin, so we can later
    # easily filter ngrams by the state index without including
    # padded elements.
    if pad_symbol:
        seq = chain((pad_symbol,)* pre_order, seq, (pad_symbol,) * post_order)
        seq = list(seq)

    # We generate the collection of ngrams for counting occurences
    # with Python itertool's `zip_longest()` function, so we can
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
    # The `zip_longest()` function will then make an iterator
    # aggregating elements from each of these iterables, using
    # `None` as a standard fill element for missing values
    # (which we can easily filter with `all()`).
    subseqs = [
        list(subseq) for subseq in
        zip_longest(*[seq[i:] for i in range(order)])
        if all(subseq)]

    print(seq)
    print(subseqs)

    # We can now collect all the skipping sequences, caching the
    # various indexes for quicker extraction. Given the way
    # Python indexes lists, we need to perform a conditional
    # extraction based inthe value of post_order (otherwise
    # we would get a `postctx_idx` of -1+1, i.e. zero, which
    # would return the entire subsequence; work-arounds would
    # make the code too complex for something so simple, and
    # this has the advantage of being faster in cases of
    # only subsequent ngrams, which are more common).
    if not post_order:
        ngrams = (
            (tuple(subseq[:-1] + [_ELEMENT]), subseq[-1], state_idx + pre_order)
            for state_idx, subseq
            in enumerate(subseqs))
    else:
        elem_idx = -1 - post_order
        prectx_idx = elem_idx - pre_order
        postctx_idx = -post_order
        ngrams = (
            # pre-context + element + post_context
            (tuple(subseq[prectx_idx:elem_idx] + [_ELEMENT] + subseq[postctx_idx:]),
             # transition
             subseq[elem_idx],
             # state index, discounting the start boundaries
             state_idx + pre_order)
            # loop
            for state_idx, subseq in enumerate(subseqs)
        )
        
    # Only keep the elements that refer to actually observed states, i.e.,
    # whose index is not None. The boundary information will be carried by
    # the states.
#    ngrams = [ngram for ngram in ngrams if ngram[2] is not None]

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
