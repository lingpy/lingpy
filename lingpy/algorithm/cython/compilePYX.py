
import os


def pyx2py(infile):

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
    newlines = []
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
        if old_line != line:
            print(old_line, line)

    with open('_'+infile+'.py','w') as f:
        f.write(''.join(newlines))
    print('converter file {0}'.format(infile))

                    


def main():

    for inf in ['calign','cluster','misc','malign','talign']:
        pyx2py(inf)
    input('pause')
    os.system('cython calign.pyx')
    os.system('cython cluster.pyx')
    os.system('cython misc.pyx')
    os.system('cython malign.pyx')
    os.system('cython talign.pyx')


    os.system('gcc -c -fPIC -I/usr/include/python3.4m/ -I/usr/lib/python3.4/site-packages/numpy/core/include/ calign.c')
    os.system('gcc -c -fPIC -I/usr/include/python3.4m/ -I/usr/lib/python3.4/site-packages/numpy/core/include/ cluster.c')
    os.system('gcc -c -fPIC -I/usr/include/python3.4m/ -I/usr/lib/python3.4/site-packages/numpy/core/include/ misc.c')
    os.system('gcc -c -fPIC -I/usr/include/python3.4m/ -I/usr/lib/python3.4/site-packages/numpy/core/include/ talign.c')
    os.system('gcc -c -fPIC -I/usr/include/python3.4m/ -I/usr/lib/python3.4/site-packages/numpy/core/include/ malign.c')

    os.system('gcc -shared calign.o -o calign.so')
    os.system('gcc -shared calign.o -o misc.so')
    os.system('gcc -shared calign.o -o malign.so')
    os.system('gcc -shared calign.o -o talign.so')
    os.system('gcc -shared calign.o -o cluster.so')

    print('Successful')


if __name__ == '__main__':
    main()

