from lingpy import *

lex = LexStat('GER.csv')

lex._get_scorer()
lex.cluster(method='lexstat',threshold=0.725)

from lingpy.evaluate.acd import diff
diff(lex)
