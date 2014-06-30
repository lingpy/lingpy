from numpy import *
from lingpy import *

lex = LexStat('../data/eurasian/samoyedic.csv')

print("Loaded " + str(len(lex.concept)) + " concepts in " + str(len(lex.language)) + " languages: " + ", ".join(lex.language));

lex.get_scorer()
