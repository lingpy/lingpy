"""
Script handles compilation of Cython files to C and also to C-Extension modules.
"""
import os
from sys import argv

def pyx2py(infile, debug=False):

    with open(infile+'.pyx') as f:
        data = f.readlines()
    
    words = [
        'bint',
        'int',
        'float',
        'list',
        'str',
        'object',
        'tuple',
        'dict',
        'double',
        ]
    newlines = ['from __future__ import unicode_literals\n']
    for line in data:
        old_line = line
        if 'cdef extern from' in line:
            line = '# [autouncomment] '+line
            line += 'from numpy import sqrt\n'
        elif 'double sqrt( double x)' in line:
            line = '# [autouncomment] ' + line
        elif not '=' in line:
            if 'cdef' in line:
                for w in words:
                    if 'cdef '+w+' ' in line:
                        line = '# [autouncomment] '+line
                        break
            else:
                for w in words:
                    if ' '+w+' ' in line:
                        line = line.replace(w+' ','')
                        break
        else:
            if 'cdef' in line:
                for w in words:
                    if 'cdef '+w+' ' in line:
                        line = line.replace('cdef '+w+' ','')
                        break
            else:
                for w in words:
                    if ' '+w+' ' in line:
                        line = line.replace(w+' ','')
                        break
        newlines += [line]
        if old_line != line and debug:
            print(old_line, line)

    with open('_'+infile+'.py','w') as f:
        f.write(''.join(newlines))

def main():
    
    scripts = ['calign','cluster','misc','malign','talign']
    for script in scripts:
        print('[i] compiling {0}...'.format(script))
        pyx2py(script)
        os.system('cython '+script+'.pyx')
        cmd = 'gcc -c -fPIC -I/usr/include/python{1}m/ {0}.c'.format(script,
                argv[1])
        os.system(cmd)
        os.system('gcc -shared '+script+'.o -o '+script+'.so')
        print('... done.')

if __name__ == '__main__':
    main()

