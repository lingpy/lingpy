"""
MPA -- Multiple Phonetic Alignments using the SCA approach.
"""
__author__ = "Johann-Mattis List"
__date__ = "2012-12-11"

# external modules
import sys
import os
from sys import argv

# lingpy-modules
from lingpy.basic.wordlist import Wordlist
from lingpy.align.multiple import Multiple

# check for args
args = {}

for arg in argv:
    try:
        name,value = arg.split('=')
        args[name] = value
    except:
        if arg.endswith('.csv'):
            args['file'] = arg

if not args or '--help' in args or '-h' in args:
    print("Usage: python[3] mpa.py file=FILENAME folder=FOLDER")
    sys.exit()
elif not 'file' in args:
    print("[i] Input file missing. Exiting...")
elif not 'folder' in args or not args['folder']:
    print("[i] Output folder missing. Exiting...")

try:
    os.mkdir(args['folder'])
except:
    pass

wl = Wordlist(args['file'])

errors = []
for c in wl.concept:

    # get all identifiers in a flat list
    idfs = wl.get_list(c,flat=True)
    
    # get all taxa
    taxa = [wl[i,'doculect'] for i in idfs]
    
    # get all strings
    strings = [wl[i,'ortho_parse'] for i in idfs]

    # get all tokens and join them for input in aligner
    tokens = [' '.join(string).replace("ˈ",'') for string in strings]

    # check tokens for double entries:
    check = 0
    for i,t in enumerate(tokens):
        if ' #' in t:
            tokens[i] = t.split(' #')[0]
            check += 1
    
    error = False
    # check for missing entries
    for i,t in enumerate(tokens):
        if not t:
            print("[!] Error in {0}...".format(idfs[i]))
            errors.append(idfs[i])
            error = True
    
    if not error:
        try:

            # align the stuff
            # print number for error-tracing
            print("[i] Aligning words for concept {0} (removed {1} double entries)...".format(c,check))
            
            # create multiple object
            msa = Multiple(tokens)

            # make the alignment using default settings
            msa.prog_align()

            # write alignment to file, cannot be done automatically at the moment, so
            # we need this workaround
            out = open(args['folder']+'/'+c.replace('"','')+'.msa','w')
            out.write('PAD\n')
            out.write(c+'\n')
            for i,taxon in enumerate(taxa):
                out.write('{0:15}\t'.format(taxon)+'\t'.join(msa.alm_matrix[i])+'\n')
            out.close()
            print("... done.")
        except:
            print("[!] Error in concept {0}".format(c))
            errors.append(c)

    


