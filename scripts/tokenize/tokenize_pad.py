# -*- coding: utf-8 -*-
"""
Test script for the orthography / tokenization module.
"""

__author__ = "Steven Moran"
__date__ = "2010-11-06"

from lingpy import *

# rcParams["debug"]=True

wl = Wordlist("pad_data_qlc.qlc")
wl.tokenize('pad_orthography_profile')
wl.output("qlc", filename="tokenized-pad_data")
