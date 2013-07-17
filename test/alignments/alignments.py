from lingpyd import *

alm = Alignments('KSL.csv')
print("[i] Testing progressive alignments")
alm.align(method='progressive')
print("[i] Testing library alignments.")
alm.align(method='library')
input()
# switch to verbose
rc(verbose=True)
alm.align(method='library',iterate=True,swap_check=True)
input()
# check pairwise module
p = Pairwise("toxta","tugatera")
p.align(distance=True)
print(p)
input()
# check multiple module
m = Multiple(["toxta","tugatera","dota"])
m.prog_align()
print(m)
input()
m.lib_align()
print(m)
input()
