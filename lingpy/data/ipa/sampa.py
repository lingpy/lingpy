# *-* coding: utf-8 *-*
# These lines were automatically added by the 3to2-conversion.
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
import types

from ...settings import *
from lingpy import log

try:
    import regex as re
    if not isinstance(re, types.ModuleType):
        # this is the case when creating the docs!
        raise ImportError
except ImportError:
    import re
    log.missing_module('regex')

import sys
import os
import codecs

# data for sampa2ipa (Peter Kleiwegs implementation)
path = os.path.split(    
        os.path.abspath(        
            __file__       
            )        
        )[0]

f = codecs.open(path + '/sampa.csv','r','utf-8')
xsdata = []
_xsKeys = [' ']
xs = {' ': ' '}

for line in f:
    line = line.strip('\n').strip('\r')
    if line and not line.startswith('#'):
        key,val = line.split('\t')
        try:
            assert key not in xs
        except:
            sys.stderr.write(key + '\n')
            sys.stderr.flush()
        _xsKeys.append(key)
        xs[key] = eval('"""'+val+'"""')

_kk = []
for _k in _xsKeys:
    _kk.append(re.escape(_k))
_kk.sort(reverse = True)  # long before short
_xsPat = '|'.join(_kk)
reXS = re.compile('(' + _xsPat + ')|(.)')
"""
The regular expression used in the sampa2unicode-converter 
is taken from an
algorithm for the conversion of XSAMPA to IPA (Unicode) by Peter 
Kleiweg <http://www.let.rug.nl/~kleiweg/L04/devel/python/xsampa.html>.
@author: Peter Kleiweg
@date: 2007/07/19
"""
