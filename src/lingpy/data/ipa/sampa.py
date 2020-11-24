# *-* coding: utf-8 *-*
"""
The regular expression used in the sampa2unicode-converter  is taken from an
algorithm for the conversion of XSAMPA to IPA (Unicode) by Peter
Kleiweg <http://www.let.rug.nl/~kleiweg/L04/devel/python/xsampa.html>.
@author: Peter Kleiweg
@date: 2007/07/19
"""
from __future__ import print_function, division, unicode_literals
import re
import sys
import codecs

from lingpy.util import data_path


# data for sampa2ipa (Peter Kleiwegs implementation)
xsdata = []
_xsKeys = [' ']
xs = {' ': ' '}

for line in codecs.open(data_path('ipa', 'sampa.csv'), 'r', 'utf-8'):
    line = line.strip('\n').strip('\r')
    if line and not line.startswith('#'):
        key,val = line.split('\t')
        if key in xs and xs[key] != val:
            raise ValueError("Keys encode too many values.")
        _xsKeys.append(key)
        xs[key] = eval('"""' + val + '"""')

_kk = []
for _k in _xsKeys:
    _kk.append(re.escape(_k))
_kk.sort(reverse=True)  # long before short
_xsPat = '|'.join(_kk)
reXS = re.compile('(' + _xsPat + ')|(.)')
