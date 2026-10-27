"""
Calculation of various diversity functions.
"""

import collections
import math


def simpson(cognates):
    """
    Calculate Simpson's diversity for a list of cognates.


    Note
    ----

    For details on the Simpson index, see the Wikipedia entry
    https://de.wikipedia.org/wiki/Simpson-Index. This entry also links to the
    original paper.
    """
    count = 0
    types = set(cognates)
    base = len(cognates) * (len(cognates) - 1)
    for cogid in types:
        prop = cognates.count(cogid)
        count += (prop * (prop - 1)) / base
    return count

    
def shannon(cognates, normalized=False):
    """
    Calculate Shannong diversity for a list of cognates.

    Parameters
    ----------
    normalized : bool

    Note
    ----
    For details on the Shannong index, see the Wikipedia entry
    https://de.wikipedia.org/wiki/Shannon-Index. This entry also links to the
    original paper. The normalization just divides the index by the maximum
    index, which is equivalent to the natural logarithm of the number of
    different cognate sets.
    """

    count = 0
    types = set(cognates)
    base = len(cognates)
    for cogid in types:
        prop = cognates.count(cogid) / base
        count += prop * math.log(prop)
    if normalized:
        return - count / math.log(base)
    return - count or 0.0
