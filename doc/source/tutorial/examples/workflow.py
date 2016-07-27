# author   : Johann-Mattis List, Steven Moran
# email    : mattis.list@gmail.com
# created  : 2013-04-11 14:34
# modified : 2013-11-12 13:27
"""
An example workflow of LingPy.

Usage: python workflow.py
"""

__author__="Johann-Mattis List, Steven Moran"
__date__="2013-11-12"

from lingpy import *

# create a lexstat object from your data
print("[i] Loading wordlist as a LexStat object.")
lex = LexStat('DOGON_tokens.qlc')

# calculate lex-stat-scorer, using default parameters
print('[i] Calculating the scorer.')
lex.get_scorer()

# calculate cognates using default parameters and the lexstat method
print('[i] Clustering the words into cognate sets.')
lex.cluster(method='lexstat',threshold=0.5,verbose=False)

# write data to file, specify the columns that shell be output
print('[i] Writing parts of the data to file.')
lex.output(
    'qlc',
    filename='DOGON_lexstat',
    subset=True,
    prettify=False,
    ignore="all",
    cols=['concepts','taxa','counterpart','tokens','lexstatid']
    )

print('[i] Loading data into the Alignments class.')
msa = Alignments('DOGON_lexstat.qlc',ref='lexstatid')

print('[i] Carrying out alignment analyses.')
msa.align(
        method='library',
        output = True        
        )

print('[i] Writing results to wordlist file.')
msa.output('qlc',filename='DOGON_alignments')

print("[i] Plotting the alignments.")
msa.output('html',filename='DOGON')

print("[i] Starting with the borrowing detection.")
from lingpy.compare.phylogeny import PhyBo
tre = PhyBo('DOGON_lexstat.qlc',
        ref='lexstatid', degree=180)

tre.analyze(runs='weighted', mode='mixed')

print("[i] Plotting the results.")
# make nice labels for the taxon-names
labels = dict(
        [(t,t[:9]) for t in tre.taxa]
        )
tre.get_MLN(tre.best_model)
tre.plot_MLN(
        tre.best_model,
        filename='DOGON',
        fileformat='svg',
        figsize=(15,7),
        labels=labels
        )
