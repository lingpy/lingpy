# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-12 13:26
# modified : 2013-07-12 13:26
"""
Module provides namespaces and data for Evolaemp applications.
"""

__author__="Johann-Mattis List"
__date__="2013-07-12"


from re import findall
from pickle import load
import os
import codecs

from ...check.exceptions import *
from ...data.model import Model,load_dvt

try:
    from ...data.derive import compile_model,compile_dvt
except ImportError:
    print(ThirdPartyModuleError('networkx').warning())

global ipa_diacritics
global ipa_vowels
global ipa_tones
global sca
global asjp
global dolgo
global art
global _color

try:
    ipa_diacritics,ipa_vowels,ipa_tones, = load_dvt(path='el')
    sca = Model('sca_el')
    asjp = Model('asjp_el')
    dolgo = Model('dolgo_el')
    art = Model('art_el')
    _color = Model('color_el')
# compile the models if they are not precompiled
except:
    compile_dvt(path='el')
    compile_model('sca_el')
    compile_model('dolgo_el')
    compile_model('art_el')
    compile_model('color_el')
    ipa_diacritics,ipa_vowels,ipa_tones = load_dvt(path='el')
    sca = Model('sca_el')
    asjp = Model('asjp_el')
    dolgo = Model('dolgo_el')
    art = Model('art_el')
    _color = Model('color_el')
