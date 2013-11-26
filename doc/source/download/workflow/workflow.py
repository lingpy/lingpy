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

# load a wordlist object
print("[i] Loading the wordlist.")
wl = Wordlist('DOGON.qlc')

# tokenize the data using the ortho-profile Heath2012.prf
print("[i] Tokenizing the wordlist.")
wl.tokenize('Heath2013.prf',column='IPA')

# output the data to the file "DOGON_tokens.csv"
print("[i] Writing wordlist to file.")
wl.output("qlc",filename="DOGON_tokens")

# create a lexstat object from the data
print("[i] Loading wordlist as a LexStat object.")
lex = LexStat('DOGON_tokens.qlc')

# calculate lex-stat-scorer, using default parameters
print('[i] Calculating the scorer.')
lex.get_scorer()

print('[i] Pickling the data.')
lex.pickle()

# calculate cognates using default parameters and the lexstat method
print('[i] Clustering the words into cognate sets.')
lex.cluster(method='lexstat',threshold=0.5,verbose=False)

# write data to file, specify the columns that shell be output
print('[i] Writing parts of the data to file.')
lex.output(
    'qlc',
    filename='DOGON_lexstat',
    subset=True,
    formatter='concepts',
    cols=['concepts','taxa','counterpart','tokens','lexstatid'],
    ignore=['scorer']
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
msa.output('html', filename='DOGON', ref='lexstatid')

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
