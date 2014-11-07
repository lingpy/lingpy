from numpy import *
from lingpy import *

lex = LexStat('samoyedic.csv')

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(lex.language)) + " languages: " + ", ".join(lex.language));

lex.get_scorer()

lex.cluster(method='lexstat', threshold=0.5)

lex.calculate('tree', ref='lexstatid', tree_calc='neighbor')

print(lex.tree.asciiArt())

lex.output('qlc', filename='samoyedic-results', subset=True, formatter='concept', cols=['concept', 'taxa', 'counterpart', 'tokens', 'lexstatid'], ignore=["scorer"])

alm = Alignments('samoyedic-results.qlc', ref='lexstatid')

alm.align(method='library', output=False)

alm.output("html", filename='samoyedic-alignments', ref='lexstatid')

print("Alignment output as HTML complete!")