from lingpyd import *

rc(verbose=True)
lex = LexStat('Witotoan.csv',check=True)
lex.get_scorer()
lex.cluster(method='lexstat')
