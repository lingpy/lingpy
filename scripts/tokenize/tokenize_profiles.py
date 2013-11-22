"""
Test script for the orthography / tokenization module.
"""

__author__ = "Steven Moran"
__date__ = "2010-11-06"

import sys
from lingpy import *

# rcParams["debug"]=True


# op and no rules - dictionary
d = Dictionary("minor1971-3-74.csv")
d.tokenize('minor1971')
# d.output("qlc", filename="tokenized-minor1971") 

# op and no rules - dictionary
d = Dictionary("burtch1983-19-262.csv")
d.tokenize('burtch1983')
# d.output("qlc", filename="tokenized-leach1969") 

# op and rules
wl = Dictionary("leach1969-67-161.csv")
wl.tokenize('leach1969')
# wl.output("qlc", filename="tokenized-leach1969") 

# op, rules
wl = Wordlist("huber1992.csv")
wl.tokenize('huber1992')
wl.output("qlc", filename="tokenized-huber1992") 

sys.exit(1)
# op, no rules
wl = Wordlist("zgraggen1980.csv")
wl.tokenize('zgraggen1980')
wl.output("qlc", filename="tokenized-zgraggen1980") 

# no op, rules
wl = Wordlist("pad_data_qlc.qlc")
wl.tokenize('pad_orthography_profile')
wl.output("qlc", filename="tokenized-pad_data")

# no op and no rules
wl = Wordlist("huber1992.csv")
wl.tokenize()
wl.output("qlc", filename="no-op-tokenization-huber1992")

# op and no rules, transform to ipa
wl = Wordlist("huber1992.csv")
wl.tokenize("huber1992-2", conversion="ipa")
wl.output("qlc", filename="ipa-tokenization-huber1992")
