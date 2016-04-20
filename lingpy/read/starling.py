"""
Basic parser for Starling data.
"""
from __future__ import unicode_literals, division, print_function

from lingpy.read.csv import csv2list
from lingpy import log
from lingpy.util import identity


def star2qlc(filename, clean_taxnames=False, debug=False):
    """
    Converts a file directly output from starling to LingPy-QLC format.
    """
    cleant = clean_taxnames or identity
    data = csv2list(filename)

    # check for strange chars in data due to notepad errors
    data[0][0] = data[0][0].replace('\ufeff', '')

    # get the header
    header = data[0]

    # debugging
    if debug:
        error = False
        log.info("Header line has length {0}.".format(len(header)))
        for line in data[1:]:
            if len(line) != len(header):  # pragma: no cover
                log.error("Error for item {0} with length {1}, expected {2}.".format(
                    '/'.join(line[0:2]), len(line), len(header)))
                error = True
        if error:  # pragma: no cover
            log.error("Errors were found, aborting function call.")
            return
        else:
            log.info("Everything went fine, carrying on with function call.")

    # determine language names in header
    taxa = []
    for i in range(len(header) - 1):
        prev = header[i]
        post = header[i + 1]

        if prev in post and '#' in post:
            taxa += [prev]

            if len(taxa) == 1:
                lngIdx = i

        if prev == 'Number':
            numIdx = i

        if prev == 'Word':
            wrdIdx = i

    log.info('starling, indices (%s, %s, %s)' % (lngIdx, numIdx, wrdIdx))
    log.info('starling, taxa: %s' % taxa)

    # start filling in the dictionary
    D = {0: [
        'DOCULECT', 'CONCEPT', 'GLOSSID', 'WORDINSOURCE', 'ORTHOGRAPHY', 'IPA', 'COGID']}

    idx = 1
    cognate_counter = 0
    current_concept = ''
    cognate_sets = []
    for line in data[2:]:
        gloss = line[wrdIdx]
        gnum = line[numIdx]

        # switch to next cognate set if there is a switch in concepts
        if current_concept != gloss and len(cognate_sets) != 0:
            max_cog = max(cognate_sets)
            cognate_counter = max_cog
            cognate_sets = []
            current_concept = gloss
        else:
            log.debug('starling, indices (%s, %s, %s)' % (
                gloss, current_concept, cognate_counter))

        for i in range(lngIdx, len(header), 2):
            word = line[i]

            if '{' in word:
                ipa = word[:word.index('{')].strip()
                ortho = word[word.index('{') + 1:word.index('}')].strip()
            else:
                ipa = word
                ortho = word

            cogid = int(line[i + 1])

            if cogid != 0 and word:
                if cogid > 0:
                    cogid = cogid + cognate_counter

                # append cognate sets, essential for raising the counter
                cognate_sets += [int(cogid)]

                taxon = cleant(header[i])

                D[idx] = [taxon, gloss, gnum, word, ortho, ipa, cogid]
                idx += 1

    # re-iterate through data and reassign cognate sets with negative ids
    for k in D:
        if k:
            cogid = D[k][-1]
            if cogid < 0:
                cogid = -cognate_counter
                cognate_counter += 1
                D[k][-1] = cogid

    return D
