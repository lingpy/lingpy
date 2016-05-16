"""
Module provides various string similarity metrics.
"""
from collections import defaultdict
import itertools as it

def ldn(a, b, normalized=True):
    """Basic Levenshtein distance without swap operation (all operations are equal costs).

    See also
    --------
    lingpy.align.pairwise.edit_dist
    lingpy.compare.strings.ldn_swap
    """
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            if a[i - 1] != b[j - 1]:
                s += 1
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la -= 1
    lb -= 1
    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def ldn_swap(a, b, normalized=True):
    """Basic Levenshtein distance with swap operation included (identifies metathesis).
    """
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            if (a[i - 1] != b[j - 1]):
                s = s + 1
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
            if i > 1 and j > 1 and a[i - 1] == b[j - 2] and a[i - 2] == b[j - 1]:
                m[i][j] = min(m[i][j], m[i - 2][j - 2] + 1)
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def bidist1(a, b, normalized=True):
    """Computes bigram-based distance.
    
    Notes
    -----
    The binary version.
    Checks if two bigrams are equal or not.
    """
    pad_symbol = "-"
    n = 2
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    m = []
    la = len(s_a) + 1
    lb = len(s_b) + 1
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            if (s_a[i - 1] != s_b[j - 1]):
                s = s + 1
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def tridist1(a, b, normalized=True):
    """Computes trigram-based distance.
    
    Notes
    -----
    The binary version.
    Checks if two trigrams are equal or not.
    """
    pad_symbol = "-"
    n = 3
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    m = []
    la = len(s_a) + 1
    lb = len(s_b) + 1
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            if (s_a[i - 1] != s_b[j - 1]):
                s = s + 1
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def bidist2(a, b, normalized=True):
    """Computes bigram based distance. 
    
    Notes
    -----
    The comprehensive version of the
    bigram distance.
    """
    pad_symbol = "-"
    n = 2

    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    m = []
    la = len(s_a)
    lb = len(s_b)
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            dis = 0.0
            dis = len([k for k in s_a[i - 1] if k not in s_b[j - 1]]) / 2.0
            # dis = dis/2.0
            if (s_a[i - 1] != s_b[j - 1]):
                s = s + dis
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def tridist2(a, b, normalized=True):
    """Computes bigram based distance. 
    
    Notes
    -----
    The comprehensive version of the
    bigram distance.
    """
    pad_symbol = "-"
    n = 3

    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]
    m = []
    la = len(s_a)
    lb = len(s_b)
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            dis = 0.0
            dis = len([k for k in s_a[i - 1] if k not in s_b[j - 1]]) / 3.0
            # dis = dis/3.0
            if (s_a[i - 1] != s_b[j - 1]):
                s = s + dis
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def bidist3(a, b, normalized=True):
    """Computes bigram based distance. 
    
    Notes
    -----
    Computes the positional version of the
    bigrams. Assigns a partial distance between two bigrams based on positional
    similarity of bigrams.
    """
    pad_symbol = "-"
    n = 2

    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]
    m = []
    la = len(s_a)
    lb = len(s_b)
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            dis = 0.0
            if s_a[i - 1][0] != s_b[j - 1][0]:
                dis += 1.0
            if s_a[i - 1][1] != s_b[j - 1][1]:
                dis += 1.0
            dis = dis / 2.0
            if (s_a[i - 1] != s_b[j - 1]):
                s = s + dis
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def tridist3(a, b, normalized=True):
    """Computes trigram based distance. 
    
    Notes
    -----
    Computes the positional version of the
    trigrams. Assigns a partial distance between two trigrams based on positional
    similarity of trigrams.
    """
    pad_symbol = "-"
    n = 3

    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]
    m = []
    la = len(s_a)
    lb = len(s_b)
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = i
    for i in range(0, lb):
        m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            dis = 0.0
            if s_a[i - 1][0] != s_b[j - 1][0]:
                dis += 1.0
            if s_a[i - 1][1] != s_b[j - 1][1]:
                dis += 1.0
            if s_a[i - 1][2] != s_b[j - 1][2]:
                dis += 1.0
            dis = dis / 3.0
            if (s_a[i - 1] != s_b[j - 1]):
                s = s + dis
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(m[la][lb])
    return float(m[la][lb]) / float(max(la, lb))


def dice(a, b, normalized=True):
    """Computes the Dice measure that measures the number of common bigrams.
    """
    la = len(a) - 1
    lb = len(b) - 1
    overlap = 0
    dicta = defaultdict(int)
    dictb = defaultdict(int)
    for i in range(len(a) - 1):
        tmp = ",".join(map(str, a[i:i + 2]))
        dicta[tmp] += 1
    for j in range(len(b) - 1):
        tmp = ",".join(map(str, b[j:j + 2]))
        dictb[tmp] += 1
    for entry in dicta:
        if (dictb.__contains__(entry)):
            overlap = overlap + min(dicta.get(entry), dictb.get(entry))
    total = la + lb
    if total == 0:
        return 0
    if not normalized:
        #        return float(2.0 * overlap)
        return float(total) - float(2.0 * overlap)
    return 1.0 - (float(2.0 * overlap) / float(total))


def lcs(a, b, normalized=True):
    """Computes the longest common subsequence between two strings.
    """
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            if (a[i - 1] == b[j - 1]):
                m[i][j] = m[i - 1][j - 1] + 1
            else:
                m[i][j] = max(m[i][j - 1], m[i - 1][j])
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def bisim1(a, b, normalized=True):
    """computes the binary version of bigram similarity.
    """
    pad_symbol = "-"
    n = 2
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            if (s_a[i - 1] == s_b[j - 1]):
                m[i][j] = m[i - 1][j - 1] + 1
            else:
                m[i][j] = max(m[i][j - 1], m[i - 1][j])
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def trisim1(a, b, normalized=True):
    """Computes the binary version of trigram similarity.
    """
    pad_symbol = "-"
    n = 3
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            if (s_a[i - 1] == s_b[j - 1]):
                m[i][j] = m[i - 1][j - 1] + 1
            else:
                m[i][j] = max(m[i][j - 1], m[i - 1][j])
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def bisim2(a, b, normalized=True):
    """Computes bigram similarity "the comprehensive version".
    
    Notes
    -----
    Computes the
    number of common 1-grams between two n-grams.
    """
    pad_symbol = "-"
    n = 2
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]
    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            sim = len([k for k in s_a[i - 1] if k in s_b[j - 1]]) / 2.0
            m[i][j] = max(m[i][j - 1], m[i - 1][j], m[i - 1][j - 1] + sim)
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def trisim2(a, b, normalized=True):
    """Computes tri-sim "the comprehensive version". 
    
    Notes
    -----
    Simply computes the number of common 1-grams between two n-grams instead of
    calling LCS as should be done in :evobib:`Kondrak2005` paper. Note that the
    LCS for a trigram can be computed in O(n) time if we asssume that list
    lookup is in constant time.
    """
    pad_symbol = "-"
    n = 3
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            sim = len([k for k in s_a[i - 1] if k in s_b[j - 1]]) / 3.0
            m[i][j] = max(m[i][j - 1], m[i - 1][j], m[i - 1][j - 1] + sim)
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def bisim3(a, b, normalized=True):
    """Computes bi-sim the positional version.
    
    Notes
    -----
    The partial similarity between two
    bigrams is defined as the number of matching 1-grams at each position.
    """
    pad_symbol = "-"
    n = 2
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            sim = 0.0
            if s_a[i - 1][0] == s_b[j - 1][0]:
                sim += 1.0
            if s_a[i - 1][1] == s_b[j - 1][1]:
                sim += 1.0
            sim = sim / 2.0
            m[i][j] = max(m[i][j - 1], m[i - 1][j], m[i - 1][j - 1] + sim)
    la = la - 1
    lb = lb - 1
    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def trisim3(a, b, normalized=True):
    """Computes tri-sim the "positional version". 
    
    Notes
    -----
    Simply computes the
    number of matching 1-grams in each position.
    """
    pad_symbol = "-"
    n = 3
    m = []
    la = len(a) + 1
    lb = len(b) + 1
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    count = max(0, len(s_a) - n + 1)
    s_a = [tuple(s_a[i:i + n]) for i in range(count)]
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)
    count = max(0, len(s_b) - n + 1)
    s_b = [tuple(s_b[i:i + n]) for i in range(count)]

    for i in range(0, la):
        m.append([])
        for j in range(0, lb):
            m[i].append(0)
        m[i][0] = 0
    for i in range(0, lb):
        m[0][i] = 0
    for i in range(1, la):
        for j in range(1, lb):
            sim = 0.0
            if s_a[i - 1][0] == s_b[j - 1][0]:
                sim += 1.0
            if s_a[i - 1][1] == s_b[j - 1][1]:
                sim += 1.0
            if s_a[i - 1][2] == s_b[j - 1][2]:
                sim += 1.0
            sim = sim / 3.0
            m[i][j] = max(m[i][j - 1], m[i - 1][j], m[i - 1][j - 1] + sim)
    la = la - 1
    lb = lb - 1

    if not normalized:
        return float(max(la, lb)) - float(m[la][lb])
    return 1.0 - (float(m[la][lb]) / float(max(la, lb)))


def jcd(a, b, normalized=True):
    """
    Computes the bigram-based Jaccard Index.
    """
    la = len(a) - 1
    lb = len(b) - 1
    overlap = 0
    dicta = defaultdict(int)
    dictb = defaultdict(int)
    for i in range(len(a) - 1):
        tmp = ",".join(map(str, a[i:i + 2]))
        dicta[tmp] += 1
    for j in range(len(b) - 1):
        tmp = ",".join(map(str, b[j:j + 2]))
        dictb[tmp] += 1
    for entry in dicta:
        if dictb.__contains__(entry):
            overlap = overlap + min(dicta.get(entry), dictb.get(entry))
    total = la + lb - overlap
    if total == 0:
        return 0
    if not normalized:
        return float(total) - float(overlap)
    return 1.0 - (float(overlap) / float(total))


def jcdn(a, b, normalized=True):
    """
    Computes the bigram and trigram-based Jaccard Index
    """
    la = len(a) - 1
    lb = len(b) - 1
    overlap = 0
    n = 3
    pad_symbol = "-"
    s_a = it.chain((pad_symbol,) * (n - 1), a)
    s_a = it.chain(s_a, (pad_symbol,) * (n - 1))
    s_a = list(s_a)
    s_b = it.chain((pad_symbol,) * (n - 1), b)
    s_b = it.chain(s_b, (pad_symbol,) * (n - 1))
    s_b = list(s_b)

    dicta = defaultdict(int)
    dictb = defaultdict(int)
    for i in range(len(s_a) - 1):
        for k in range(1, n + 1):
            tmp = ",".join(map(str, s_a[i:i + k]))
            dicta[tmp] += 1
    for j in range(len(s_b) - 1):
        for k in range(1, n + 1):
            tmp = ",".join(map(str, s_b[i:i + k]))
            dictb[tmp] += 1
    for entry in dicta:
        if dictb.__contains__(entry):
            overlap = overlap + min(dicta.get(entry), dictb.get(entry))
    total = la + lb - overlap
    if total == 0:
        return 0
    if not normalized:
        return float(total) - float(overlap)
    return 1.0 - (float(overlap) / float(total))


def prefix(a, b, normalized=True):
    """Computes the longest common prefix between two strings.
    """
    la = len(a)
    lb = len(b)
    minl = min(la, lb)
    maxl = max(la, lb)
    pref = 0
    for i in range(minl):
        if a[i] == b[i]:
            pref += 1
    if not normalized:
        return float(maxl) - float(pref)
    return 1.0 - (float(pref) / float(maxl))


def xdice(a, b, normalized=True):
    """Computes the skip 1 character version of Dice.
    """
    la = len(a) - 2
    lb = len(b) - 2
    overlap = 0
    dicta = defaultdict(int)
    dictb = defaultdict(int)
    for i in range(len(a) - 2):
        tmp = ",".join(map(str, [a[i], a[i + 2]]))
        dicta[tmp] += 1
    for j in range(len(b) - 2):
        tmp = ",".join(map(str, [b[j], b[j + 2]]))
        dictb[tmp] += 1
    for entry in dicta:
        if dictb.__contains__(entry):
            overlap = overlap + min(dicta.get(entry), dictb.get(entry))
    total = la + lb

    if total == 0 or la < 1 or lb < 1:
        return 0
    if not normalized:
        return float(total) - float(2.0 * overlap)
    return 1.0 - (float(2 * overlap) / float(total))


def trigram(a, b, normalized=True):
    """Computes the number of common trigrams between two strings.
    """
    la = len(a) - 2
    lb = len(b) - 2
    overlap = 0
    dicta = defaultdict(int)
    dictb = defaultdict(int)
    for i in range(len(a) - 2):
        tmp = ",".join(map(str, a[i:i + 3]))
        dicta[tmp] += 1
    for j in range(len(b) - 2):
        tmp = ",".join(map(str, b[j:j + 3]))
        dictb[tmp] += 1
    for entry in dicta:
        if (dictb.__contains__(entry)):
            overlap = overlap + min(dicta.get(entry), dictb.get(entry))
    total = la + lb

    if total == 0 or la < 1 or lb < 1:
        return float(1.0)
    if not normalized:
        return float(total) - float(2.0 * overlap)
    return 1.0 - (float(2 * overlap) / float(total))


def ident(a, b):
    """Computes the identity between two strings. If yes, returns 1, else, returns 0.
    """
    return 1 if a.__eq__(b) else 0

def xxdice(a, b, normalized=True):
    """Returns the XXDice between two strings. 
    
    Notes
    -----
    Taken from :evobib:`Brew1996`.

    """
    la = len(a) - 1
    lb = len(b) - 1
    overlap = 0
    dicta = defaultdict(list)
    dictb = defaultdict(list)
    for i in range(len(a) - 1):
        tmp = ",".join(map(str, a[i:i + 2]))
        dicta[tmp].append(i)
    for j in range(len(b) - 1):
        tmp = ",".join(map(str, b[j:j + 2]))
        dictb[tmp].append(j)
    for entry in dicta:
        if (dictb.__contains__(entry)):
            pos_a = dicta[entry]
            pos_b = dictb[entry]
            for m, n in it.product(pos_a, pos_b):
                overlap += 1.0 / (1.0 + (m - n) ** 2)
    total = la + lb
    if total == 0:
        return 0
    if not normalized:
        return float(total) - float(2.0 * overlap)
    return 1.0 - (float(2.0 * overlap) / float(total))
