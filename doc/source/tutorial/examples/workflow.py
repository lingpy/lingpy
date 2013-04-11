# author   : Johann-Mattis List, Steven Moran
# email    : mattis.list@gmail.com
# created  : 2013-04-11 14:34
# modified : 2013-04-11 14:34
"""
An example workflow of LingPy.

Usage: python workflow.py
"""

__author__="Johann-Mattis List, Steven Moran"
__date__="2013-04-11"


from lingpy import *

# load a wordlist object
print("[i] Loading the wordlist.")
wl = Wordlist('DOGON.csv')

# tokenize the data using the ortho-profile Heath2012.prf
print("[i] Tokenizing the wordlist.")
wl.tokenize('Heath2012.prf')

# output the data to the file "DOGON_tokens.csv"
print("[i] Writing wordlist to file.")
wl.output("csv",filename="DOGON_tokens")

# create a lexstat object from the data
print("[i] Loading wordlist as a LexStat object.")
lex = LexStat('DOGON_tokens.csv')

# calculate lex-stat-scorer, using default parameters
print('[i] Calculating the scorer.')
lex.get_scorer(verbose=False)

print('[i] Pickling the data.')
lex.pickle()

# calculate cognates using default parameters and the lexstat method
print('[i] Clustering the words into cognate sets.')
lex.cluster(method='lexstat',threshold=0.5,verbose=False)

# write data to file, specify the columns that shell be output
print('[i] Writing parts of the data to file.')
lex.output(
    'csv',
    filename='DOGON_lexstat',
    subset=True,
    formatter='concepts',
    cols=['concepts','taxa','counterpart','tokens','lexstatid']
    )

print('[i] Loading data into the Alignments class.')
msa = Alignments('DOGON_lexstat.csv',cognates='lexstatid')

print('[i] Carrying out alignment analyses.')
msa.align(
        method='library',
        output = True,
        verbose = False
        )

print('[i] Writing results to wordlist file.')
msa.output('csv',filename='DOGON_alignments')
msa.output('alm',filename='DOGON',cognates='lexstatid')

print("[i] Plotting the alignments.")
from lingpy.convert.plot import alm2html
alm2html('DOGON.alm',filename='DOGON')

print("[i] Carrying out the analysis for borrowing detection.")
from lingpy.convert.plot import alm2html
alm2html('DOGON.alm',filename='DOGON')

from lingpy.compare.borrowing.trebor import TreBor
tre = TreBor('DOGON_lexstat.csv',cognates='lexstatid')

tre.analyze(runs='weighted',mode='mixed')

print("[i] Plotting the results.")
tre.plot_MLN(filename='DOGON',fileformat='svg')
