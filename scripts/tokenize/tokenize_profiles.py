"""
Test script for the orthography / tokenization module.
"""

__author__ = "Steven Moran"
__date__ = "2010-11-06"

from lingpy import *

rcParams["debug"]=True

# op, rules
wl = Wordlist("huber1992.csv")
wl.tokenize('huber1992')
wl.output("qlc", filename="tokenized-huber1992") 

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
