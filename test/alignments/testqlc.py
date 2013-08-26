from lingpyd import *

from lingpyd.basic.wordlist import WordlistTest

wl = WordlistTest('KSL.csv')
wl.calculate('tree')
wl.pickle()
rc(verbose=True)
wl = WordlistTest('KSL.csv')

wl.output('qlc',filename='wordlisttest')
