
import os

os.system('cython calignx.pyx')
input('ok')
os.system('gcc -c -fPIC -I/usr/include/python3.3m/ -I/usr/lib/python3.3/site-packages/numpy/core/include/ calignx.c')
os.system('gcc -shared calignx.o -o calignx.so')

print('Successful')
