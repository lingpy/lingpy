#! /usr/bin/env python
from __future__ import division,print_function
from lingpy.data.derive import compile_model
from scipy.spatial.distance import squareform
from time import sleep
from pickle import dump

asjp = {}

score = open('score','r').read()
score = score.split('\n')
del score[-1]

dicto = {}
for line in score:
    lin = line.split('\t')
    dicto[lin[0]] = lin[1:]

letters = []
for i in range(len(score)):
    score[i] = score[i].split('\t')
    letters.append(score[i][0])
    del score[i][0]

matrix = []
for i in range(len(score)):
    for l in letters:
        if i < len(dicto[l]):
            matrix.append(float(dicto[l][i]))

matrix = squareform(matrix)
consonants = ['p'] + letters
consonant_matrix = matrix.copy()

score = open('vows_score','r').read()
score = score.split('\n')
del score[-1]

dicto = {}
for line in score:
    lin = line.split('\t')
    dicto[lin[0]] = lin[1:]

letters = []
for i in range(len(score)):
    score[i] = score[i].split('\t')
    letters.append(score[i][0])
    del score[i][0]

matrix = []
for i in range(len(score)):
    for l in letters:
        if i < len(dicto[l]):
            matrix.append(float(dicto[l][i]))

matrix = squareform(matrix)

vowel_matrix = matrix.copy()
vowels = ['i'] + letters

for i in range(len(vowel_matrix)):
    vowel_matrix[i][i] = 40

for i in range(len(consonant_matrix)):
    consonant_matrix[i][i] = 40

for i in range(31):
    for j in range(31):
        asjp[consonants[i],consonants[j]] = consonant_matrix[i][j]

for i in range(7):
    for j in range(7):
        asjp[vowels[i],vowels[j]] = vowel_matrix[i][j]

for l in vowels:
    asjp[l,'X'] = 0
    asjp['X',l] = 0

for l in consonants:
    asjp[l,'X'] = 0
    asjp['X',l] = 0

asjp['X','X'] = 0
for v in vowels:
    for c in consonants:
        asjp[v,c] = -20
        asjp[c,v] = -20

for key in asjp.keys():
    if asjp[key] == 0:
        asjp[key] = 0
    else:
        asjp[key] = int(asjp[key]+0.5)

for v1 in vowels:
    for v2 in vowels:
        asjp[v1,v2] = int(asjp[v1,v2] * 0.25 + 0.5) + 10
asjp['i','y'] = -2
asjp['y','i'] = -2
asjp['u','w'] = -2
asjp['w','u'] = -2
asjp['u','v'] = -4
asjp['v','u'] = -4
asjp['u','f'] = -6
asjp['f','u'] = -6

keys = []
for keyA,keyB in asjp.keys():
    keys.append((keyA,keyB))

for keyA,keyB in keys:
    asjp[keyA,'+'] = -20
    asjp['+',keyB] = -20
    asjp[keyA,'0'] = 0
    asjp['0',keyB] = 0


asjp['X','+'] = -5
asjp['+','X'] = -5
asjp['+','+'] = 0 # swaps
asjp['0','0'] = 0 # missing values
asjp['X','0'] = 0
asjp['0','X'] = 0

for i in '0123456':
    for j in '0123456':
        if i == j:
            asjp[i,j] = 10
        else:
            asjp[i,j] = 5

keys = []
for keyA,keyB in asjp.keys():
    keys.append((keyA,keyB))


for keyA,keyB in keys:
    for i in '123456':
        if keyA not in '123456' and keyB not in '123456':
            asjp[keyA,i] = -20
            asjp[i,keyB] = -20
    asjp[keyA,'_'] = -50
    asjp['_',keyB] = -50

asjp['_','_'] = 0

for x in asjp.keys():
    asjp[x] = asjp[x] / 4.0
    if asjp[x] > 0 and asjp[x] != 10:
        asjp[x] += 0.75 * asjp[x]
    elif asjp[x] < 0:
        asjp[x] += 0.75 * asjp[x]

out = open('scorer.bin','wb')
dump(asjp,out)
out.close()
compile_model('asjp')
print("[i] Compilation of the ASJP model was successful!")
sleep(1)
