#! /usr/bin/env python
# *-* coding:utf-8 *-*

from pickle import dump

converter = {}

alphabet = 'abcdefghijklmnopqrstuvwxyz'
alphabet += alphabet.upper()
alphabet += ' '

for x in alphabet:
    converter[x] = x

scorer = {}
for x in alphabet:
    for y in alphabet:
        if x == y and 'X' not in (x,y):
            scorer[x,y] = 1
        elif x != y and 'X' not in (x,y):
            scorer[x,y] = -1
        else:
            scorer[x,y] = -1

scorer['X','X'] = 0

out = open('converter.bin','wb')
dump(converter,out)
out.close()
out = open('scorer.bin','wb')
dump(scorer,out)
out.close()
raw_input('finished...')
