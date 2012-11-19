import sys
import os
import regex as re

# data for sampa2ipa (Peter Kleiwegs implementation)
path = os.path.split(    
        os.path.abspath(        
            __file__       
            )        
        )[0]

f = open(path + '/sampa.csv')
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
