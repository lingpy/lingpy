# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-04 13:58
# modified : 2013-03-07 18:02
"""
This script tests various SCA routines.
"""

__author__="Johann-Mattis List"
__date__="2013-03-07"

from lingpy import *

# define a set of qlc-strings
qlc_strings = [
    "# v o l d e m o r t #",
    "# w a l d e m a r #",
    "# v l a d i m i r #",
    ]

# create an msa-object
msa = Multiple(qlc_strings)

# carry out an SCA alignment, progressive style
msa.prog_align()

# print results to terminal
print("")
print("Progressive Alignment, SCA-Algorithm")
print(msa)

# print results to terminal
print("")
print("Progressive Alignment, Basic-Algorithm")
msa.prog_align(classes=False,sonar=False,gop=-1)
print(msa)

# check for swapped sites
msa.prog_align()
swap = msa.swap_check(swap_penalty=-1)
if swap:
    print("")
    print("Detected swap of columns {0[0][0]} and {0[0][2]}".format(msa.swap_index))

# check the iterations
#try:
msa.iterate_orphans()
msa.iterate_clusters(0.5)
msa.iterate_similar_gap_sites()
msa.iterate_all_sequences()
print("")
print("Iteration-test was successful.")
#except:
#    print("Error in iteration test.")

# carry out pairwise tests
msa.get_pairwise_alignments()
pairs = [(''.join(a[0]).replace('-',''),''.join(a[1]).replace('-','')) for a in msa.alignments.values()]
# carrying out pairwise alignments
print("")
print("Pairwise Alignments")
pw = Pairwise(pairs)
pw.align(pprint=True,distance=True)


# load sca-file
msa = SCA('data/test.msq')
msa.lib_align()
print("")
print("Multiple alignment with MSA from text-file.")
print("Writing result to file.")
msa.output('msa',filename='data/test_out')

print("")
print("Testing pairwise alignment using SCA.")
psa = SCA('data/test.psq')
psa.align(pprint=True,distance=True)

