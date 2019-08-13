from __future__ import division
from lingpy.sequence.sound_classes import tokens2class, prosodic_string, tokens2morphemes
from lingpy.align.multiple import mult_align
from collections import defaultdict
from itertools import combinations

def _scorer():
    scoredict = {}
    for c1, c2 in combinations(list('CcVvT'), r=2):
        scoredict[c1, c2] = 0
        scoredict[c2, c1] = 0
    for c in 'CcVvT':
        scoredict[c, c] = 5
    return scoredict

def pattern_consensus(patterns, scoredict):

    alms = mult_align(patterns, 
            scoredict=scoredict)
    cons = []
    for j in range(len(alms[0])):
        col = [alms[i][j] for i in range(len(alms))]
        sounds = sorted(
                set([x for x in col if x != '-']), 
                key=lambda x: col.count(x)) 
        token = '|'.join(sounds)
        if '-' in col:
            cons += ['('+token+')']
        else:
            cons += [token]
    return cons

def cv_templates(
        wordlist, language, segments='tokens', converter=None, cutoff=0.1,
        output='markdown', examples=3, scoredict=None, splitter=False
        ):
    """Create CV templates from wordlist data."""
    templates = defaultdict(list)
    idxs = wordlist.get_list(col=language, flat=True)
    sounds = defaultdict(list)

    def str_(list_):
        return ', '.join([' '.join(l) for l in list_[:examples]])
    
    if not converter:
        converter = lambda x: prosodic_string(x, _output='CcV') 
    scoredict = scoredict or _scorer()
    if not splitter:
        splitter = lambda x: filter(None, tokens2morphemes(x))

    for idx in idxs:
        segs = wordlist[idx, segments]
        for word in splitter(segs):
            cv = converter(word)
            templates[cv] += [word]
            for sound, symbol in zip(word, cv):
                 sounds[sound, symbol] += [word]

    # retrieve percentile 
    lengths = sum([len(v) for v in templates.values()])
    perc = lengths - (cutoff * lengths)
    
    patterns, ignored = [], []
    score = 0
    for k, v in sorted(templates.items(), key=lambda x: len(x[1]), reverse=True):
        l = len(v)
        if score + l > perc:
            ignored += [[k, l, v]]
        else:
            patterns += [[k, l, v]]
        score += l
    
    # compute pattern consensus
    consensus = pattern_consensus([list(p[0]) for p in patterns], scoredict)

    # extract initials
    sound_table = []
    for k, v in sorted(sounds.items(), key=lambda x: (x[0][1], len(x[1]))):
        sound_table += [(k[0], k[1], len(v), v)]
        
    if output == 'markdown':
        out = 'Pattern | Frequency | Examples\n --- | --- | --- \n'
        score = 0
        for i, (p, l, v) in enumerate(patterns):
            out += '{0:15} | {1:5} | {2}\n'.format(p, l, str_(v))
            score += l
        count = 1
        out += '\nSound | Context | Frequency | Examples\n --- | --- | --- | --- \n'
        for sound, context, l, vals in sound_table:
            out += '{0} | {1} | {2} | {3} \n'.format(
                sound, context, l, str_(vals))
            
        out += '\n* **coverage:** {0} out of {1} patterns in the data\n'.format(
                score, lengths)
        out += '* **pattern consensus:** {0}\n'.format(' '.join(consensus))
        return out

    return patterns, ignored, sound_table

    

