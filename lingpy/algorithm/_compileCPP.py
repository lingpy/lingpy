#! /usr/bin/env python

import os

#a = raw_input('Pfad: ')
#b = raw_input('Name: ')

#if a:
#    path = a + '/' + b
#else:
#    path = b
path = 'align'
cmd1 = 'g++ -c -fPIC -I/usr/include/python3.3m ' + path + '.cpp'
cmd2 = 'g++ -shared ' + path + '.o -o '+ path + '.so -lpython3.3m -lboost_python3'
print cmd1
os.system(cmd1)
print cmd2
os.system(cmd2)
raw_input("Press enter to quit.")

