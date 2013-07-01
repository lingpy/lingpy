# -*- coding: utf-8 -*-
#!/usr/bin/env python3


import sys, os
import zipfile
import csv
import codecs
import glob

import lingpy

def main(argv):

    if len(argv) < 2:
        print("call: extract_componentas_wordlist.py component")
        sys.exit(1)

    component = argv[1]

    if not os.path.exists("sources.csv"):
        try:
            import requests
        except:
            print("Module 'requests' not found. I cannot download the data "
                  "automatically for you.\n\nPlease download manually at:\n"
                  "http://www.quanthistling.info/data/downloads/csv/data.zip")
            sys.exit(1)

        r = requests.get(
            "http://www.quanthistling.info/data/downloads/csv/data.zip")
        with open("data.zip", "wb") as f:
            f.write(r.content)

        z = zipfile.ZipFile("barwiki.zip")
        z.extractall()

    sources = csv.reader(codecs.open("sources.csv", "r", "utf-8"), delimiter="\t")
    witotoan_sources = list()
    for source in sources:
        if source[5] == component and source[1] == "dictionary": # and source[3] == "True"
            witotoan_sources.append(source[0])

    concepts = lingpy.spanish_swadesh_list()
    cm = lingpy.ConceptComparerSpanishStem()
    cg = lingpy.ConceptGraph(concepts, "spa", cm)

    for f in glob.glob("*.csv"):
        if ("-" in f and f[:f.index("-")] in witotoan_sources) or \
                ("." in f and f[:f.index(".")] in witotoan_sources):
            print("Adding {0}...".format(f))
            di = lingpy.Dictionary(f)
            cg.add_dictionary(di)

    cg.output_wordlist("{0}.csv".format(component))

if __name__ == "__main__":
    main(sys.argv)