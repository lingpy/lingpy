# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-12 13:26
# modified : 2013-07-12 13:26
"""
Module provides namespaces and data for Evolaemp applications.
"""

__author__="Johann-Mattis List"
__date__="2013-07-12"

from ...settings import rcParams,rc
rc(schema='evolaemp')

ipa_diacritics = rcParams['diacritics']
ipa_vowels = rcParams['vowels']
ipa_tones = rcParams['tones']
sca = rcParams['sca']
asjp = rcParams['asjp']
dolgo = rcParams['dolgo']
art = rcParams['art']
_color = rcParams['_color']
