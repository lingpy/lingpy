"""
Package provides generic workflow modules for LingPy.
"""
from __future__ import division, print_function, unicode_literals
import json
from functools import partial

from lingpy.settings import rcParams
from lingpy.align.sca import Alignments
from lingpy.align._align.confidence import get_correspondences
from lingpy.compare.lexstat import LexStat
from lingpy import util


class Workflow(object):
    """
    Class provides access to generic workflows.

    Parameters
    ----------
    infile : str
        A tsv-file providing the input data for the given workflow.
    """

    def __init__(self, infile):

        # we don't do anything specific here, we just assign the input file as
        # an attribute of the Workflow class
        self.infile = infile

    def cognate_detection(self, **keywords):
        """
        Method runs a cognate detection analysis.
        """
        kw = dict(
            align_method='progressive',
            align_mode=rcParams['align_mode'],
            align_modes=rcParams['align_modes'],
            cluster_method=rcParams['lexstat_cluster_method'],
            cognate_method='sca',
            cognate_mode='overlap',
            defaults=False,
            factor=rcParams['align_factor'],
            gap_weight=rcParams['gap_weight'],
            gop=rcParams['align_gop'],
            iteration=False,
            lexstat_modes=rcParams['lexstat_modes'],
            limit=rcParams['lexstat_limit'],
            merge_vowels=rcParams['merge_vowels'],
            model=rcParams['sca'],
            export="html",
            preprocessing=False,
            preprocessing_method=rcParams['lexstat_preprocessing_method'],
            preprocessing_threshold=rcParams['lexstat_preprocessing_threshold'],
            rands=rcParams['lexstat_rands'],
            ratio=rcParams['lexstat_ratio'],
            ref="customid",
            restricted_chars=rcParams['restricted_chars'],
            restriction='',
            runs=rcParams['lexstat_runs'],
            scale=rcParams['align_scale'],
            scoring_method=rcParams['lexstat_scoring_method'],
            swap_check=False,
            threshold=rcParams['lexstat_threshold'],
            tree_calc=rcParams['align_tree_calc'],
            vscale=rcParams['lexstat_vscale'],
            outfile=False,
            sonar=True,
        )

        # first load
        kw.update(keywords)
        if kw['defaults']:
            return kw

        # carry out lexstat cluster analysis
        self.lex = LexStat(self.infile, **kw)

        # reset filename if it is not defined
        kw['outfile'] = kw['outfile'] or self.lex.filename + '_lingpy'

        # check for traditional lexstat analysis
        if kw['cognate_method'] == 'lexstat':
            self.lex.get_scorer(
                method=kw['scoring_method'], modes=kw['lexstat_modes'], **kw)

        self.lex.cluster(method=kw['cognate_method'], mode=kw['cognate_mode'], **kw)

        # align the data
        self.alms = Alignments(self.lex, **kw)
        kw['scoredict'] = self.lex.cscorer \
            if kw['cognate_method'] == 'lexstat' else self.lex.bscorer

        self.alms.align(
            method=kw['align_method'],
            mode=kw['align_mode'],
            modes=kw['align_modes'],
            **kw)

        if 'tsv' in kw['export']:
            self.alms.output(
                'tsv',
                filename=kw['outfile'],
                ignore=['scorer', 'json', 'taxa', 'msa'],
                **kw)
        if 'html' in kw['export']:
            corrs, occs = get_correspondences(self.alms, kw['ref'])

            # serialize the wordlist
            wl = {}
            for concept in self.alms.concepts:
                entries = self.alms.get_list(concept=concept, flat=True)
                cogids = [self.alms[idx, kw['ref']] for idx in entries]
                words = [self.alms[idx, 'ipa'] for idx in entries]
                alms = [self.alms[idx, 'alignment'] for idx in entries]
                langs = [self.alms[idx, 'doculect'] for idx in entries]

                checkalm = lambda x: x if type(x) == str else ' '.join(x)

                wl[concept] = [list(k) for k in sorted(
                    zip(
                        langs,
                        [str(x) for x in entries],
                        words,
                        [str(x) for x in cogids],
                        [checkalm(x) for x in alms],
                    ),
                    key=lambda x: int(x[3]))]

                # make simple gloss id for internal use as id
                gloss2id = list(
                    zip(
                        self.alms.concepts,
                        [str(x) for x in range(1, len(self.alms.concepts) + 1)]))
                id2gloss = dict([[b, a] for a, b in gloss2id])
                gloss2id = dict(gloss2id)

                txt = ''
                txt += 'CORRS = ' + json.dumps(corrs) + ';\n'
                txt += 'LANGS = ' + json.dumps(self.alms.taxa) + ';\n'
                txt += 'OCCS = ' + json.dumps(occs) + ';\n'
                txt += 'WLS = ' + json.dumps(wl) + ';\n'
                txt += 'GlossId = ' + json.dumps(gloss2id) + ';\n'
                txt += 'IdGloss = ' + json.dumps(id2gloss) + ';\n'
                txt += 'FILE = "' + kw['outfile'] + '.tsv";\n'

                tpath = partial(util.data_path, 'templates')

                tname = 'jcov.{0}.html'.format(
                    'remote' if 'remote' in kw['export'] else 'direct')
                content = util.read_text_file(tpath(tname))

                util.write_text_file(
                    kw['outfile'] + '.html',
                    content.format(
                        CORRS=txt,
                        JCOV=util.read_text_file(tpath('jcov.js')),
                        STYLE=util.read_text_file(tpath('jcov.css')),
                        VENDOR=util.read_text_file(tpath('jcov.vendor.js')),
                        DIGHL=util.read_text_file(tpath('jcov.dighl.js'))))
