# -*- coding: utf-8 -*-
#!/usr/bin/env python3


import sys, os
import zipfile
import csv
import codecs
import glob

import lingpy.meaning.concepts

def main(argv):

    if len(argv) < 2:
        print("call: extract_componentas_wordlist2.py component [use_profiles]")
        sys.exit(1)

    component = argv[1]
    use_profiles = False
    if len(argv) > 2 and argv[2] == "use_profiles":
        use_profiles = True

    # use non-stemmed Spanish swadesh list
    concepts = concepts = lingpy.meaning.concepts.spanish_swadesh_list(False)
    output_file = lingpy.meaning.concepts.extract_component_as_wordlist(
        component, use_profiles, concepts, "spa")

    print("Written wordlist {0}".format(output_file))

if __name__ == "__main__":
    main(sys.argv)