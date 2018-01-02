# *-* coding: utf-8 *-*
"""
Basic functions for the conversion from LingPy to CLDF and vice versa.
"""
from __future__ import unicode_literals
import unicodedata
from collections import defaultdict
from clldutils.path import Path, read_text
from csvw.metadata import TableGroup
from clldutils.dsv import reader, Dialect
from clldutils.misc import slug

from lingpy import util
from lingpy.basic.wordlist import Wordlist
from lingpy.convert.html import template_path

try:
    from pycldf import Wordlist as CLDF_Wordlist
    cldf = True
except ImportError:
    cldf = False


def from_cldf(path, to=Wordlist):
    """
    Load data from CLDF into a LingPy Wordlist object or similar.

    Parameters
    ----------
    path : str
        The path to the metadata-file of your CLDF dataset.
    to : ~lingpy.basic.wordlist.Wordlist
        A ~lingpy.basic.wordlist.Wordlist object or one of the descendants
        (LexStat, Alignmnent).

    Note
    ----
    This function does not offer absolute flexibility regarding the data you
    can input so far. However, it can regularly read CLDF-formatted data into
    LingPy and thus allow you to use CLDF data in LingPy analyses.

    Todo
    ----
    Add support for partial cognates.
    """
    tbg = TableGroup.from_file(path)
    forms = tbg.tabledict['forms.csv']

    # obtain the dictionaries to convert ids to values
    taxa = {t['ID']: (t['Name'], t['glottocode']) for t in
            tbg.tabledict['languages.csv']}
    concepts = {c['ID']: (c['Name'], c['Concepticon_ID']) for c in 
            tbg.tabledict['parameters.csv']}

    # create dictionary
    D = {}
    id2idx = {}
    for i, row in enumerate(forms):
        # check for numeric ID
        if row['ID'].isdigit:
            idx = int(row['ID'])
        else:
            idx = i+1
        id2idx[row['ID']] = idx

        doculect, glottocode = taxa[row['Language_ID']]
        concept, concepticon_id = concepts[row['Parameter_ID']]
        D[idx] = [doculect, glottocode, concept, concepticon_id] + [row.get(f,
            '') or '' for f in ['form_in_source', 'Form', 'Segments',
                'Comment', 'Source']]
    # add the header
    D[0] = ['doculect', 'glottocode', 'concept', 'concepticon_id', 'value',
            'form', 'tokens', 'note', 'source']

    # convert to wordlist (simplifies handling)
    wordlist = to(D)
    
    # add cognates if they are needed
    if 'cognates.csv' in tbg.tabledict:
        cognates = {id2idx[row['Form_ID']]: (row['Cognateset_ID'],
            row['Alignment']) for row in tbg.tabledict['cognates.csv']}
        wordlist.add_entries('cogid', cognates, lambda x: x[0] or 0)
        wordlist.add_entries('alignment', cognates, lambda x: x[1] or '')

    return wordlist


def to_cldf(wordlist, path='cldf', source_path=None, ref="cogid",
        segments="tokens", form="ipa", note='note', form_in_source="value",
        source=None, alignment=None):
    """Convert a wordlist in LingPy to CLDF.
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        A regular Wordlist object (or similar).
    path : str (default='cldf')
        The name of the directory to which the files will be written.
    source_path : str (default=None)
        If available, specify the path of your BibTex file with the sources.
    ref : str (default="cogid")
        The column in which the cognate sets are stored.
    segments : str (default="tokens")
        The column in which the segmented phonetic strings are stored.
    form : str (default="ipa")
        The column in which the unsegmented phonetic strings are stored.
    note : str (default=None)
        The column in which you store your comments.
    form_in_source : str (default=None)
        The column in which you store the original form in the source.
    source : str (default=None)
        The column in which you store your source information. 
    alignment : str (default="alignment")
        The column in which you store the alignments.
    """
    if not cldf:
        raise ValueError('The package pycldf needs to be installed')

    # create cldf-dataset
    ds = CLDF_Wordlist.in_dir(path)
    # add sources if they are available
    ds.add_sources(
            read_text(source_path) if source_path else '')
    # add components
    ds.add_component('LanguageTable')
    ds.add_component('ParameterTable', 'Concepticon_ID')
    ds.add_component('CognateTable')
    ds.add_columns('FormTable', 'form_in_source')

    languages, parameters, forms, cognates = {}, {}, [], []
    for idx in wordlist:
        lid = slug(wordlist[idx, 'doculect'])
        if lid not in languages:
            languages[lid] = dict(
                    ID=lid,
                    Name=wordlist[idx, 'doculect'],
                    glottocode = wordlist[idx, 'glottocode'])

        pid = wordlist[idx, 'concepticon_id'] or slug(wordlist[idx, 'concept'])
        if pid not in parameters:
            parameters[pid] = dict(
                ID=pid,
                Name=wordlist[idx, 'concept'],
                Concepticon_ID=wordlist[idx, 'concepticon_id'])

        forms.append(dict(
            ID=str(idx),
            Language_ID=lid,
            Parameter_ID=pid,
            form_in_source=wordlist[idx, form_in_source] or '' if form_in_source else '',
            Form=wordlist[idx, form] or '' if form else '',
            Segments=wordlist[idx, segments] or '' if segments else '',
            Source=[wordlist[idx, source]] or [] if source else [],
            Comment=wordlist[idx, note] or '' if note else ''))

        if ref:
            cognates.append(dict(ID=str(idx), Form_ID=str(idx),
                Cognateset_ID=wordlist[idx, ref], Alignment=wordlist[idx,
                    alignment] or [''] if alignment else ['']))

    ds.write(
        FormTable=forms,
        LanguageTable=languages.values(),
        ParameterTable=parameters.values(),
        CognateTable=cognates)


