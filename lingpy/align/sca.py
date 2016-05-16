# *-* coding: utf-8 *-*
"""
Basic module for pairwise and multiple sequence comparison.

The module consists of four classes which deal with pairwise and multiple
sequence comparison from the *sequence* and the *alignment* perspective. The
sequence perspective deals with unaligned sequences. The *alignment*
perspective deals with aligned sequences.

"""
from __future__ import print_function, division, unicode_literals
import os
from collections import Counter, defaultdict

from six import text_type

from lingpy import __version__
from lingpy.read.qlc import read_msa, normalize_alignment, reduce_alignment
from lingpy.settings import rcParams
from lingpy.basic.wordlist import Wordlist
from lingpy.convert import html
from lingpy.convert.strings import msa2str
from lingpy.sequence.sound_classes import (
    ipa2tokens, tokens2class, class2tokens, prosodic_string, prosodic_weights,
    tokens2morphemes,
)
from lingpy.align.multiple import Multiple
from lingpy.align.pairwise import Pairwise
from lingpy.algorithm import misc
from lingpy.align._align import confidence
from lingpy import util
from lingpy import log


class MSA(Multiple):
    """
    Basic class for carrying out multiple sequence alignment analyses.

    Parameters
    ----------
    infile : file
        A file in ``msq``-format or ``msa``-format.
    merge_vowels : bool (default=True)
        Indicate, whether neighboring vowels should be merged into
        diphtongs, or whether they should be kept separated during the
        analysis.
    comment : char (default='#')
        The comment character which, inserted in the beginning of a line,
        prevents that line from being read.
    normalize : bool (default=True)
        Normalize the alignment, that is, add gap characters for all sequences
        which are shorter than the longest sequence, and delete all columns
        from the alignment in which only gaps occur.


    Examples
    --------
    Get the path to a file from the testset.

    >>> from lingpy import *
    >>> path = rc("test_path")+'harry.msq'

    Load the file into the Multiple class.

    >>> mult = Multiple(path)

    Carry out a progressive alignment analysis of the sequences.

    >>> mult.prog_align()

    Print the result to the screen:

    >>> print(mult)
    w    o    l    -    d    e    m    o    r    t
    w    a    l    -    d    e    m    a    r    -
    v    -    l    a    d    i    m    i    r    -

    Notes
    -----
    There are two possible input formats for this class: the MSQ-format, and
    the MSA-format (see :ref:`msa_formats` for details). This class directly
    inherits all methods of the :py:class:`~lingpy.align.multiple.Multiple`
    class.

    """

    def __init__(self, infile, **keywords):
        util.setdefaults(
            keywords,
            comment=rcParams['comment'],
            diacritics=rcParams['diacritics'],
            vowels=rcParams['vowels'],
            tones=rcParams['tones'],
            combiners=rcParams['combiners'],
            breaks=rcParams['breaks'],
            stress=rcParams['stress'],
            merge_vowels=rcParams['merge_vowels'],
            ids=False,
            header=True,
            normalize=True)
        self.comment = keywords['comment']

        # initialization checks first, whether we are dealing with msa-files or
        # with other, unaligned, sequence files and starts the
        # loading-procedures accordingly
        if isinstance(infile, dict):
            self._init_dict(infile, **keywords)
        else:
            if infile.endswith('.msa') or infile.endswith('.msq'):
                tmp = read_msa(infile, **keywords)
                self._init_dict(tmp, **keywords)

    def _init_dict(self, initdict, **keywords):
        """
        Initialize by passing a dictionary with the relevant values.
        """
        for key in initdict:
            if key not in 'seqs':
                setattr(self, key, initdict[key])

        Multiple.__init__(self, initdict['seqs'], **keywords)

        if 'alignment' in initdict:
            self.alm_matrix = initdict['alignment']
        if 'local' in initdict:
            self.local = initdict['local']
        if 'swaps' in initdict:
            self.swaps = initdict['swaps']

    def ipa2cls(self, **keywords):
        """
        Retrieve sound-class strings from aligned IPA sequences.

        Parameters
        ----------
        model : str (default='sca')
            The sound-class model according to which the sequences shall be
            converted.

        Notes
        -----
        This function is only useful when an ``msa``-file with already
        conducted alignment analyses was loaded.
        """
        util.setdefaults(keywords, model=rcParams['sca'], stress=rcParams['stress'])
        self.classes = []
        self.model = keywords['model']

        # redefine the sequences of the Multiple class
        class_strings = [
            tokens2class(seq.split(' '), self.model, stress=keywords['stress'])
            for seq in self.seqs]

        # define the scoring dictionaries according to the methods
        aligned_seqs = [alm for alm in self.alm_matrix]
        for i in range(len(aligned_seqs)):
            self.classes.append(list(
                ''.join(
                    class2tokens(class_strings[i], aligned_seqs[i])).replace('-', 'X')))

    def output(
            self,
            fileformat='msa',
            filename=None,
            sorted_seqs=False,
            unique_seqs=False,
            **keywords):
        """
        Write data to file.

        Parameters
        ----------
        fileformat : { "psa", "msa", "msq" }
            Indicate which data should be written to file. Select between:

            * "psa" -- output of all pairwise alignments in ``psa``-format,
            * "msa" -- output of the multiple alignment in ``msa``-format, or
            * "msq" -- output of the multiple sequences in ``msq``-format.
            * "html" -- output of the multiple alignment in ``html``-format.

        filename : str
            Select a specific name for the outfile, otherwise, the name of
            the infile will be taken by default.

        sorted_seqs : bool
            Indicate whether the sequences should be sorted or not (applys only
            to 'msa' and 'msq' output.

        unique_seqs : bool
            Indicate whether only unique sequences should be written to file or
            not.

        """
        util.setdefaults(keywords, wordlist=False, timestamp=False)

        if fileformat in ['html', 'tex']:
            with util.TemporaryPath(suffix='.msa') as tmp:
                self.output(
                    fileformat='msa',
                    filename=os.path.splitext(tmp)[0],
                    sorted_seqs=sorted_seqs,
                    unique_seqs=unique_seqs)
                if 'filename' not in keywords:
                    keywords['input_file'] = os.path.split(self.infile)[1]
                    keywords['filename'] = filename

                getattr(html, 'msa2' + fileformat)(tmp, **keywords)
                return

        # create a specific format string in order to receive taxa of equal length
        mtax = max([len(t) for t in self.taxa])
        txf = '{0:.<' + text_type(mtax) + '}'

        with util.TextFile((filename or self.infile) + '.' + fileformat) as out:
            # start writing data to file
            out.write(self.dataset + '\n')

            if fileformat in ['msq', 'msa']:
                out.write(self.seq_id + '\n')

            if not sorted_seqs or fileformat == 'psa':
                for i, taxon in enumerate(self.taxa):
                    if fileformat == 'msq':
                        out.write(txf.format(taxon) + '\t' + self.seqs[i] + '\n')
                    elif fileformat == 'msa':
                        out.write(txf.format(taxon) + '\t')
                        out.write('\t'.join(self.alm_matrix[i]) + '\n')
                    elif fileformat == 'psa':
                        if not hasattr(self, 'alignments'):
                            self.get_pairwise_alignments(new_calc=False)
                        else:
                            pass
                        for j, taxonB in enumerate(self.taxa):
                            if i < j:
                                try:
                                    almA, almB, score = self.alignments[i, j]
                                except:
                                    almB, almA, score = self.alignments[j, i]

                                out.write('{0} ({1}, {2})\n'.format(
                                    self.seq_id,
                                    taxon,
                                    taxonB))
                                out.write(txf.format(taxon) + '\t')
                                out.write('\t'.join(almA) + '\n')
                                out.write(txf.format(taxonB) + '\t')
                                out.write('\t'.join(almB) + '\n')
                                out.write('{0} {1:.2f}\n\n'.format(self.comment, score))
            elif sorted_seqs:
                if fileformat == 'msa':
                    alms = ['\t'.join(alm) for alm in self.alm_matrix]
                else:
                    alms = [seq for seq in self.seqs]

                if not unique_seqs:
                    taxalms = zip(self.taxa, alms)
                    taxalms = sorted(taxalms, key=lambda x: x[1])
                elif unique_seqs:
                    uniqs = sorted([x[0] for x in self.uniseqs.values()])
                    taxa = [self.taxa[x] for x in uniqs]
                    alms = [alms[x] for x in uniqs]
                    taxalms = zip(taxa, alms)
                    taxalms = sorted(taxalms, key=lambda x: x[1])

                for taxon, alm in taxalms:
                    if fileformat not in ['html', 'tex']:
                        out.write(txf.format(taxon) + '\t' + alm + '\n')

            if fileformat == 'msa':
                if hasattr(self, 'local'):
                    if self.local:
                        out.write(txf.format("LOCAL") + '\t')
                        tmp = ['.'] * len(self.alm_matrix[0])
                        for i in self.local:
                            tmp[i] = '*'
                        out.write('\t'.join(tmp) + '\n')

                if hasattr(self, 'swaps'):
                    if self.swaps:
                        out.write(txf.format('SWAPS') + '\t')
                        tmp = ['.'] * len(self.alm_matrix[0])
                        for i in self.swaps:
                            tmp[i[0]] = '+'
                            tmp[i[1]] = '-'
                            tmp[i[2]] = '+'
                        out.write('\t'.join(tmp) + '\n')
                if hasattr(self, 'merge'):
                    if len(set(self.merge.values())) < len(set(self.merge.keys())):
                        out.write(txf.format('MERGE') + '\t')
                        tmp = ['.'] * len(self.alm_matrix[0])
                        start = False
                        before = 1
                        for k in sorted(self.merge):
                            if self.merge[k] == before and not start:
                                tmp[k - 1] = '<'
                                start = True
                            if self.merge[k] == before and start:
                                tmp[k] = '-'
                            if self.merge[k] != before and start:
                                start = False
                                tmp[k - 1] = '>'
                                before += 1
                            elif self.merge[k] != before and not start:
                                before += 1
                        out.write('\t'.join(tmp) + '\n')
                if hasattr(self, 'proto'):
                    out.write(txf.format('PROTO') + '\t')
                    out.write('\t'.join(self.proto) + '\n')
                if hasattr(self, 'consensus'):
                    out.write(txf.format("CONSE") + '\t')
                    out.write('\t'.join(self.consensus) + '\n')

            if keywords['timestamp']:
                out.write('# Created using LingPy\n')
                if hasattr(self, 'params'):
                    out.write('# Parameters: ' + self.params + '\n')
                out.write('# Created: {0}\n'.format(rcParams['timestamp']))


class PSA(Pairwise):
    """
    Basic class for dealing with the pairwise alignment of sequences.

    Parameters
    ----------
    infile : file
        A file in ``psq``-format.
    merge_vowels : bool (default=True)
        Indicate, whether neighboring vowels should be merged into
        diphtongs, or whether they should be kept separated during the
        analysis.
    comment : char (default='#')
        The comment character which, inserted in the beginning of a line,
        prevents that line from being read.

    Attributes
    ----------
    taxa : list
        A list of tuples containing the taxa of all sequence pairs.
    seqs : list
        A list of tuples containing all sequence pairs.
    tokens : list
        A list of tuples containing all sequence pairs in a tokenized form.

    Notes
    -----
    In order to read in data from text files, two different file formats can be
    used along with this class: the PSQ-format, and the PSA-format (see
    :ref:`psa_formats` for details). This class inherits the methods of the
    :py:class:`~lingpy.align.pairwise.Pairwise` class.

    """
    def __init__(self, infile, **keywords):
        util.setdefaults(
            keywords,
            comment=rcParams['comment'],
            diacritics=rcParams['diacritics'],
            vowels=rcParams['vowels'],
            tones=rcParams['tones'],
            combiners=rcParams['combiners'],
            breaks=rcParams['breaks'],
            stress=rcParams["stress"],
            merge_vowels=rcParams['merge_vowels'],
        )
        # add comment-char
        self.comment = keywords['comment']

        self.infile, suffix = os.path.splitext(os.path.basename(infile))

        # import the data from the input file
        data = []
        if not suffix and os.path.exists(infile + '.psq'):
            infile = infile + '.psq'

        for line in util.read_text_file(infile, lines=True):
            if not line.startswith(self.comment):
                data.append(line)

        # set the first parameters
        # delete the first line of the data, since they are no longer needed
        self.dataset = data.pop(0)

        # append the other lines of the data, they consist of triplets,
        # separated by double line breaks
        self.taxa = []
        self.pairs = []
        self.seq_ids = []

        # check the ending of the infile
        if suffix == '.psa':
            self.alignments = []
            handle_data = self._handle_psa_data
        else:
            handle_data = self._handle_seq_data

        i = 0
        while i <= len(data) - 3:
            try:
                self.seq_ids.append(data[i])
                handle_data(data, i)
                i += 4
            except:
                log.warn("Line {0} of the data is probably miscoded.".format(i + 1))
                i += 1

        self.pair_num = len(self.pairs)
        Pairwise.__init__(self, self.pairs, **keywords)

    def _handle_psa_data(self, data, i):
        """
        Load a ``psa``-file.
        """

        almA = data[i + 1].split('\t')
        almB = data[i + 2].split('\t')
        taxonA = almA.pop(0)
        taxonB = almB.pop(0)

        kw = dict(condition=lambda k: k != '-')
        self.taxa.append((taxonA, taxonB))
        self.pairs.append((util.dotjoin(*almA, **kw), util.dotjoin(*almB, **kw)))
        self.alignments.append(
            ([text_type(a) for a in almA], [text_type(b) for b in almB], 0))

    def _handle_seq_data(self, data, i):
        """
        Load a ``psq``-file.
        """
        taxonA, seqA = data[i + 1].split('\t')
        taxonB, seqB = data[i + 2].split('\t')
        self.taxa.append((taxonA.strip('.'), taxonB.strip('.')))
        self.pairs.append((seqA, seqB))

    def output(self, fileformat='psa', filename=None, **keywords):
        """
        Write the results of the analyses to a text file.

        Parameters
        ----------
        fileformat : { 'psa', 'psq' }
            Indicate which data should be written to file. Select between:

            * 'psa' -- output of all pairwise alignments in ``psa``-format,
            * 'psq' -- output of the multiple sequences in ``psq``-format.

        filename : str
            Select a specific name for the outfile, otherwise, the name of
            the infile will be taken by default.

        """
        assert fileformat in ['psa', 'psq']
        util.setdefaults(
            keywords,
            gop=-2,
            model=rcParams['sca'],
            transform=rcParams['align_transform'],
            scores=False)
        filename = filename or self.infile

        # define the outfile and check, whether it already exists
        outfile = filename + '.' + fileformat
        # check whether outfile already exists
        if os.path.isfile(outfile):
            outfile = filename + '_out.' + fileformat

        with util.TextFile(outfile) as out:
            out.write(self.dataset + '\n')
            # if data is simple, just write simple data to file
            if fileformat == 'psq':
                for i, (a, b) in enumerate(self.pairs):
                    # determine longest taxon in order to create a format string
                    # for taxa of equal length
                    txf = '{0:.<' + text_type(max([len(t) for t in self.taxa[i]])) + '}'

                    out.write(self.seq_ids[i] + '\n')
                    out.write(txf.format(self.taxa[i][0]) + '\t' + a + '\n')
                    out.write(txf.format(self.taxa[i][1]) + '\t' + b + '\n\n')
            else:  # if fileformat == 'psa':
                for i, (a, b, c) in enumerate(self.alignments):
                    # determine longest taxon in order to create a format string
                    # for taxa of equal length
                    txf = '{0:.<' + text_type(max([len(t) for t in self.taxa[i]])) + '}'

                    out.write(self.seq_ids[i] + '\n')
                    out.write(txf.format(self.taxa[i][0]) + '\t' + '\t'.join(a) + '\n')
                    out.write(txf.format(self.taxa[i][1]) + '\t' + '\t'.join(b) + '\n')

                    if keywords['scores']:
                        # get partial alignment scores
                        scores = []
                        idxA, idxB = 0, 0
                        proA = self.weights[i][0]
                        proB = self.weights[i][1]

                        for x, y in zip(a, b):
                            if '-' not in (x, y):
                                scores += [self.model(x, y)]
                                idxA += 1
                                idxB += 1
                            else:
                                if x == '-':
                                    scores += [keywords['gop'] * proB[idxB]]
                                    idxB += 1
                                elif y == '-':
                                    scores += [keywords['gop'] * proA[idxA]]
                                    idxA += 1

                        out.write(txf.format(self.comment) + '\t' + '\t'.join(
                            ['{0:.2f}'.format(s) for s in scores]
                        ) + '\t{0:.2f}\n'.format(sum(scores)))
                    out.write('{0} {1:.2f}'.format(self.comment, c) + '\n\n')


class Alignments(Wordlist):
    """
    Class handles Wordlists for the purpose of alignment analyses.

    Parameters
    ----------
    infile : str
        The name of the input file that should conform to the basic format of
        the `~lingpy.basic.wordlist.Wordlist` class and define a specific ID
        for cognate sets.
    row : str (default = "concept")
        A string indicating the name of the row that shall be taken as the
        basis for the tabular representation of the word list.
    col : str (default = "doculect")
        A string indicating the name of the column that shall be taken as the
        basis for the tabular representation of the word list.
    conf : string (default='')
        A string defining the path to the configuration file.
    ref : string (default='cogid')
        The name of the column that stores the cognate IDs.
    modify_ref : function (default=False)
        Use a function to modify the reference. If your cognate identifiers
        are numerical, for example, and negative values are assigned as
        loans, but you want to suppress this behaviour, just set this
        keyword to "abs", and all cognate IDs will be converted to their
        absolute value.

    Attributes
    ----------
    msa : dict
        A dictionary storing multiple alignments as dictionaries which can be
        directly opened and aligned with help of the ~lingpy.align.sca.SCA 
        function. The alignment objects are referenced by a key which is
        identical with the "reference" (ref-keyword) of the alignment, that is
        the name of the column which contains the cognate identifiers.


    Notes
    -----
    This class inherits from :py:class:`~lingpy.basic.wordlist.Wordlist` and
    additionally creates instances of the
    :py:class:`~lingpy.align.multiple.Multiple` class for all cognate sets that
    are specified by the *ref* keyword.

    """
    def __init__(
            self,
            infile,
            row='concept',
            col='doculect',
            conf='',
            modify_ref=False,
            _interactive=True,
            ref="cogid",
            **keywords):
        # keywords, "strings" locates, where the reference for the alignments
        # is to be found
        kw = {"segments": "tokens", "alignment": "alignment", "transcription":
                "ipa", "ref": "cogid"}
        kw.update(keywords)

        # initialize the wordlist
        Wordlist.__init__(self, infile, row, col, conf)
        self._interactive = _interactive
        self._alignment = kw['alignment'] if kw['alignment'] in \
            self.header else self._alias[kw['alignment']]
        self._segments = kw['segments'] if kw['segments'] in self.header else \
                self._alias[kw['segments']]
        self._ref = ref if ref in self.header else self._alias[ref]
        self._transcription = kw['transcription'] if kw['transcription'] in \
                self.header else self._alias[kw['transcription']]

        # check whether fuzzy (partial) alignment or normal alignment is
        # carried out, if a new namespace is used, we assume it to be plain
        self._mode = 'fuzzy' if ref in self._class_string and self._class_string[ref] \
                not in ['str', 'int'] else 'plain'
        # store loan-status
        self._modify_ref = modify_ref

        if self._segments not in self._header:
            if self._transcription in self.header:
                self.add_entries(self._segments, self._transcription,
                        ipa2tokens)
            elif self._alignment in self.header:
                self.add_entries(self._segments, self._alignment, 
                        lambda x: ' '.join([y for y in x if y not in
                            rcParams['gap_symbol']]))
            else:
                raise ValueError("No valid source for segments could be found.")
        
        self.etd = {}
        self.add_alignments(ref=self._ref, modify_ref=modify_ref)

    def add_alignments(self, ref=False, modify_ref=False):
        """
        Function adds a new set of alignments to the data.
        """
        ref = ref or self._ref
        # check for cognate-id or alignment-id in header
        try:
            self.etd[ref] = self.get_etymdict(ref=ref, modify_ref=modify_ref)
        except:
            raise ValueError("Did not find a cognate ID in the input file.")

        stridx = self.header[self._segments]
        # create the alignments by assembling the ids of all sequences
        if 'msa' not in self._meta:
            self._meta['msa'] = {ref: {}}
        if ref not in self._meta['msa']:
            self._meta['msa'][ref] = {}
        if not self._meta['msa'][ref]:
            for key, value in self.etd[ref].items():
                tmp = [x for x in value if x != 0]
                seqids = []
                for t in tmp:
                    seqids += t
                if len(seqids) > 1:
                    # set up the dictionary
                    d = {'ID': [], 'taxa': [], 'seqs': [], 'alignment': []}
                    d['dataset'] = os.path.split(os.path.splitext(self.filename)[0])[1]
                    if 'concept' in self.header:
                        concept = self[seqids[0], 'concept']
                        d['seq_id'] = '{0} ("{1}")'.format(key, concept)
                    else:
                        d['seq_id'] = '{0}'.format(key)

                    # set up the data
                    for seq in seqids:
                        if self._alignment in self.header:
                            this_string = self[seq][self.header[self._alignment]]
                        else:
                            this_string = self[seq][stridx]
                        if isinstance(this_string, text_type):
                            this_string = this_string.split(' ')
                        # check for partial cognates
                        if self._mode == 'fuzzy':
                            # split the string into morphemes
                            # FIXME add keywords for morpheme segmentation
                            morphemes = tokens2morphemes(this_string)
                            # get the position of the morpheme
                            midx = self[seq][self.header[ref]].index(key)
                            this_string = morphemes[midx]

                        d['ID'].append(seq)
                        d['taxa'].append(self[seq, 'taxa'])
                        d['seqs'].append(this_string)
                        d['alignment'].append(this_string)

                    d['alignment'] = normalize_alignment(d['alignment'])
                    self._meta['msa'][ref][key] = d

    def reduce_alignments(self, alignment=False, ref=False):
        """
        Function reduces alignments which contain columns that are marked to be \
                ignored by the user.

        Notes
        -----
        
        This function changes the data only internally: All alignments are
        checked as to whether they contain data that should be ignored. If this
        is the case, the alignments are then reduced, and stored in a specific
        item of the alignment string. If the method doesn't find any instances
        for reduction, it still makes the copies of the alignments in order to
        guarantee that the alignments with with we want to work are at the same
        place in the dictionary.
        """
        alignments = alignment or self._alignment
        ref = ref or self._ref

        if alignment not in self.header:
            raise ValueError(
                'No alignments found in your data. ' +
                'You should carry out an alignment analysis first!')

        # dictionary to add new alignments class afterwards for providing quick
        # access
        D = {}

        for k, d in self._meta['msa'][ref].items():
            ralms = reduce_alignment(d[alignment])
            if len(ralms[0]) != len(d[alignment][0]):
                log.warn('Found an alignment that could be reduced.')
            d['_'+alignment] = ralms
            for idx, alm in zip(d['ID'], d['_'+alignment]):
                D[idx] = alm
        for k in self:
            if k not in D:
                D[k] = ['']
        self.add_entries('_'+alignment, D, lambda x: x)

    def _msa2col(self, ref=False, alignment=False):
        """
        Add alignments to column (space-separated) in order to make it easy to
        parse them in the wordlist editor.
        """
        ref = ref or self._ref
        aligment = alignment or self._alignment

        tmp = {}
        # plain mode, that means, no partial alignments
        if self._mode == 'plain':
            for key, msa in self.msa[ref].items():
                for i, idx in enumerate(msa['ID']):
                    try:
                        tmp[idx] = ' '.join(msa['alignment'][i])
                    except KeyError:
                        log.error("There are no alignments in your data.  Aborting...")
                        return

            missing = [idx for idx in self if idx not in tmp]
            for m in missing:
                if self._segments in self.header:
                    tmp[m] = self[m, self._segments]
                elif self._transcription in self.header:
                    tmp[m] = ipa2tokens(self[m, self._transcription])
                elif self._alignment in self.header:
                    tmp[m] = self[m, self._alignment]
                else:
                    raise ValueError(
                        "There are no phonetic sequences (TOKENS, ALIGNMENT, or IPA) " +
                        "in your data.")
        else:
            # in this mode, we need to trace the order of the bits that make up
            # the alignemnts
            for key in self:
                # get the cognate IDs
                cogids = self[key, rcParams['ref']]
                # get the alignment
                tmp[key] = []
                for i, cogid in enumerate(cogids):
                    if cogid in self.msa[rcParams['ref']]:
                        msa = self.msa[rcParams['ref']][cogid]
                        idx = msa['ID'].index(key)
                        tmp[key] += msa['alignment'][idx]
                    else:
                        tmp[key] += tokens2morphemes(self[key, self._segments])[i]
                    # add morpheme separator as long as we don't add the last
                    # element
                    if i < len(cogids) - 1:
                        tmp[key] += [rcParams['morpheme_separator']]

        self.add_entries(alignment, tmp, lambda x: x)

    def align(self, **keywords):
        """
        Carry out a multiple alignment analysis of the data.

        Parameters
        ----------
        method : { "progressive", "library" } (default="progressive")
            Select the method to use for the analysis.
        iteration : bool (default=False)
            Set to c{True} in order to use iterative refinement methods.
        swap_check : bool (default=False)
            Set to c{True} in order to carry out a swap-check.
        model : { 'dolgo', 'sca', 'asjp' }
            A string indicating the name of the :py:class:`Model \
            <lingpy.data.model>` object that shall be used for the analysis.
            Currently, three models are supported:

            * "dolgo" -- a sound-class model based on :evobib:`Dolgopolsky1986`,

            * "sca" -- an extension of the "dolgo" sound-class model based on
              :evobib:`List2012b`, and

            * "asjp" -- an independent sound-class model which is based on the
              sound-class model of :evobib:`Brown2008` and the empirical data
              of :evobib:`Brown2011` (see the description in
              :evobib:`List2012`.

        mode : { 'global', 'dialign' }
            A string indicating which kind of alignment analysis should be
            carried out during the progressive phase. Select between:

            * "global" -- traditional global alignment analysis based on the
              Needleman-Wunsch algorithm :evobib:`Needleman1970`,

            * "dialign" -- global alignment analysis which seeks to maximize
              local similarities :evobib:`Morgenstern1996`.

        modes : list (default=[('global',-2,0.5),('local',-1,0.5)])
            Indicate the mode, the gap opening penalties (GOP), and the gap
            extension scale (GEP scale), of the pairwise alignment analyses
            which are used to create the library.

        gop : int (default=-5)
            The gap opening penalty (GOP) used in the analysis.

        scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the exension of the basic alignment
            algorithm for affine gap penalties :evobib:`Gotoh1982`.

        factor : float (default=1)
            The factor by which the initial and the descending position shall
            be modified.

        tree_calc : { 'neighbor', 'upgma' } (default='upgma')
            The cluster algorithm which shall be used for the calculation of
            the guide tree. Select between ``neighbor``, the Neighbor-Joining
            algorithm (:evobib:`Saitou1987`), and ``upgma``, the UPGMA
            algorithm (:evobib:`Sokal1958`).

        gap_weight : float (default=0)
            The factor by which gaps in aligned columns contribute to the
            calculation of the column score. When set to 0, gaps will be
            ignored in the calculation. When set to 0.5, gaps will count half
            as much as other characters.

        restricted_chars : string (default="T")
            Define which characters of the prosodic string of a sequence
            reflect its secondary structure (cf. :evobib:`List2012b`) and
            should therefore be aligned specifically. This defaults to "T",
            since this is the character that represents tones in the prosodic
            strings of sequences.
        """
        kw = dict(
            alignment=False,
            classes=rcParams['classes'],
            defaults=False,
            factor=rcParams['align_factor'],
            filename=self.filename,
            gap_weight=rcParams['gap_weight'],
            gop=rcParams['align_gop'],
            iteration=False,
            method='progressive',
            mode=rcParams['align_mode'],
            model=rcParams['sca'],
            modes=rcParams['align_modes'],
            output=False,
            plots=False,
            ref=False,
            restricted_chars=rcParams['restricted_chars'],
            scale=rcParams['align_scale'],
            scoredict=rcParams['scorer'],
            show=False,
            sonar=rcParams['sonar'],
            style='plain',
            swap_check=False,
            tree_calc=rcParams['align_tree_calc'],
        )
        kw.update(keywords)
        kw['ref'] = kw['ref'] or self._ref
        kw['alignment'] = kw['alignment'] or self._alignment

        if kw['defaults']:
            return kw

        if text_type(kw['model']) == kw['model']:
            kw['model'] = rcParams[kw['model']]

        # create a params attribute
        params = '_'.join([
            kw['method'],
            kw['model'].name,
            text_type(kw['gop']),
            '{0:.1f}'.format(kw['scale']),
            '{0:.1f}'.format(kw['factor']),
            kw['tree_calc'],
            '{0:.1f}'.format(kw['gap_weight']),
            kw['restricted_chars']
        ])

        for key, value in sorted(self.msa[kw['ref']].items(), key=lambda x: x[0]):
            log.debug("Analyzing cognate set number {0}.".format(key))

            # check for scorer keyword
            if not kw['scoredict']:
                m = SCA(value, **kw)
            else:
                # get the tokens
                numbers = [self[idx, 'numbers'] for idx in value['ID']]
                if kw['sonar']:
                    sonars = [self[idx, 'sonars'] for idx in value['ID']]
                else:
                    sonars = False
                value['seqs'] = numbers
                m = SCA(value, **kw)
                kw['sonars'] = sonars
                kw['classes'] = False

            if kw['method'] == 'progressive':
                m.prog_align(**kw)
            elif kw['method'] == 'library':
                m.lib_align(**kw)

            if kw['iteration']:
                m.iterate_similar_gap_sites()
                m.iterate_clusters(0.5)
                m.iterate_orphans()

            if kw['swap_check']:
                m.swap_check()

            # convert back to external format, if scoredict is set
            if kw['scoredict']:
                for i, alm in enumerate(m.alm_matrix):
                    tk = self[m.ID[i], 'tokens']
                    new_tk = class2tokens(tk, alm)
                    m.alm_matrix[i] = new_tk

            if hasattr(m, 'swaps'):
                self._meta['msa'][kw['ref']][key]['swaps'] = m.swaps

            self._meta['msa'][kw['ref']][key]['alignment'] = m.alm_matrix
            self._meta['msa'][kw['ref']][key]['_sonority_consensus'] = \
                m._sonority_consensus
            self._meta['msa'][kw['ref']][key]['stamp'] = rcParams['align_stamp'].format(
                m.dataset, m.seq_id, __version__, rcParams['timestamp'], params)

        self._msa2col(ref=kw['ref'], alignment=kw['alignment'])

    def get_confidence(self, scorer, ref="lexstatid", gap_weight=0.25):
        """
        Function creates confidence scores for a given set of alignments.

        Parameters
        ----------
        scorer : :py:class:`~lingpy.algorithm._misc.ScoreDict`
            A *ScoreDict* object which gives similarity scores for all segments in
            the alignment.
        ref : str (default="lexstatid")
            The reference entry-type, referring to the cognate-set to be used for
            the analysis.
        gap_weight : {loat} (default=1.0)
            Determine the weight assigned to matches containing gaps.

        """
        corrs = confidence.get_confidence(self, scorer, ref, gap_weight)
        log.info("Successfully calculated confidence values for alignments.")
        return corrs

    def __len__(self):
        return len(self.msa)

    def _plot(self, fileformat='html', **keywords):
        """
        Make an HTML plot of the aligned data.
        """
        util.setdefaults(
            keywords,
            title='LingPy - Automatic Cognate Judgments and Alignments',
            shorttitle="LingPy",
            dataset=self.filename,
            show=False,
            filename=self.filename,
            ref=False,
            confidence=False)
        keywords['ref'] = keywords['ref'] or self._ref

        with util.TemporaryPath(suffix='.alm') as tmp:
            self.output(
                'alm',
                ref=keywords['ref'],
                filename=os.path.splitext(tmp)[0],
                confidence=keywords['confidence'])
            html.alm2html(tmp, **keywords)

    def get_consensus(
            self,
            tree=False,
            gaps=False,
            classes=False,
            consensus='consensus',
            counterpart='ipa',
            weights=[],
            **keywords):
        """
        Calculate a consensus string of all MSAs in the wordlist.

        Parameters
        ----------
        msa : {c{list} ~lingpy.align.multiple.Multiple}
            Either an MSA object or an MSA matrix.
        tree : {c{str} ~lingpy.thirdparty.cogent.PhyloNode}
            A tree object or a Newick string along which the consensus shall be
            calculated.
        gaps : c{bool} (default=False)
            If set to c{True}, return the gap positions in the consensus.
        classes : c{bool} (default=False)
            Specify whether sound classes shall be used to calculate the consensus.
        model : ~lingpy.data.model.Model
            A sound class model according to which the IPA strings shall be
            converted to sound-class strings.

        """
        util.setdefaults(
            keywords, model=rcParams['sca'], gap_scale=1.0, ref=rcParams['ref'])

        # switch ref
        if keywords['ref'] != rcParams['ref']:
            rcParams['ref'] = keywords['ref']

        # reassing ref for convenience
        ref = keywords['ref']

        # check for existing alignments
        test = list(self.msa[ref].keys())[0]
        if 'alignment' not in self.msa[ref][test]:
            log.error(
                "No alignments could be found. You should carry out"
                " an alignment analysis first!")
            return

        # go on with the analysis
        cons_dict = {}
        with util.ProgressBar('CONSENSUS', len(self.etd[ref])) as progress:
            for cog in self.etd[ref]:
                progress.update()

                if cog in self.msa[ref]:
                    log.debug("Analyzing cognate set number '{0}'...".format(cog))

                    # temporary solution for sound-class integration
                    if classes == True:
                        _classes = []
                        if weights:
                            keywords['weights'] = prosodic_weights(
                                prosodic_string(self.msa[ref][cog]['_sonority_consensus'])
                            )
                        else:
                            keywords['weights'] = [
                                1.0 for i in range(len(self.msa[ref][cog]['alignment']))]

                        for alm in self.msa[ref][cog]['alignment']:
                            cls = [c for c in tokens2class(
                                alm,
                                keywords['model']
                            ) if c != '0']
                            cls = class2tokens(cls, alm)
                            _classes.append(cls)
                        _classes = misc.transpose(_classes)
                    else:
                        _classes = classes

                    cons = get_consensus(
                        self.msa[ref][cog]['alignment'],
                        classes=_classes,
                        tree=tree,
                        gaps=gaps,
                        taxa=[text_type(taxon.replace("(", "").replace(")", ""))
                              for taxon in self.msa[ref][cog]['taxa']],
                        **keywords)
                    self.msa[ref][cog]["consensus"] = cons
                # if there's no msa for a given cognate set, this set is a singleton
                else:
                    cons = self[
                        [k[0] for k in self.etd[ref][cog] if k != 0][0], counterpart]

                # add consensus to dictionary
                cons_dict[cog] = cons

        # add the entries
        self.add_entries(
            consensus, ref, lambda x: cons_dict[x], override=not self._interactive)

    def output(self, fileformat, **keywords):
        """
        Write wordlist to file.

        Parameters
        ----------
        fileformat : {"tsv", "msa", "tre", "nwk", "dst", "taxa", "starling", "paps.nex",
            "paps.csv" "html"}
            The format that is written to file. This corresponds to the file
            extension, thus 'tsv' creates a file in tsv-format, 'dst' creates
            a file in Phylip-distance format, etc. Specific output is created
            for the formats "html" and "msa":

            * "msa" will create a folder containing all alignments of all
              cognate sets in "msa"-format
            * "html" will create html-output in which words are sorted
              according to meaning, cognate set, and all cognate words are
              aligned
        filename : str
            Specify the name of the output file (defaults to a filename that
            indicates the creation date).
        subset : bool (default=False)
            If set to c{True}, return only a subset of the data. Which subset
            is specified in the keywords 'cols' and 'rows'.
        cols : list
            If *subset* is set to c{True}, specify the columns that shall be
            written to the csv-file.
        rows : dict
            If *subset* is set to c{True}, use a dictionary consisting of keys
            that specify a column and values that give a Python-statement in
            raw text, such as, e.g., "== 'hand'". The content of the specified
            column will then be checked against statement passed in the
            dictionary, and if it is evaluated to c{True}, the respective row
            will be written to file.
        ref : str
            Name of the column that contains the cognate IDs if 'starling' is
            chosen as an output format.
        missing : { str, int } (default=0)
            If 'paps.nex' or 'paps.csv' is chosen as fileformat, this character
            will be inserted as an indicator of missing data.
        tree_calc : {'neighbor', 'upgma'}
            If no tree has been calculated and 'tre' or 'nwk' is chosen as
            output format, the method that is used to calculate the tree.
        threshold : float (default=0.6)
            The threshold that is used to carry out a flat cluster analysis if
            'groups' or 'cluster' is chosen as output format.
        style : str (default="id")
            If "msa" is chosen as output format, this will write the alignments
            for each msa-file in a specific format in which the first column
            contains a direct reference to the word via its ID in the wordlist.
        ignore : { list, "all" }
            Modifies the output format in "tsv" output and allows to ignore
            certain blocks in extended "tsv", like "msa", "taxa", "json", etc.,
            which should be passed as a list. If you choose "all" as a plain
            string and not a list, this will ignore all additional blocks and
            output only plain "tsv".
        prettify : bool (default=True)
            Inserts comment characters between concepts in the "tsv" file
            output format, which makes it easier to see blocks of words
            denoting the same concept. Switching this off will output the file
            in plain "tsv".

        See also
        --------
        ~lingpy.basic.wordlist.Wordlist.output
        ~lingpy.compare.lexstat.LexStat.output

        """
        kw = dict(
            ref=rcParams['ref'],
            filename=rcParams['filename'],
            style="id",
            defaults=False,
            confidence=False)
        kw.update(keywords)
        if kw['defaults']:
            return kw

        # check for html fileformat
        if fileformat == 'html':
            self._plot(**keywords)

        # define two vars for convenience
        ref = kw['ref']
        filename = kw['filename']

        if 'cognates' in kw:
            log.deprecated('cognates', 'ref')
            ref = kw['cognates']

        if ref != rcParams['ref']:
            rcParams['ref'] = ref

        if fileformat not in ['alm', 'msa']:
            return self._output(fileformat, **kw)

        if fileformat == 'alm':
            # define the string to which the stuff is written
            out = self.filename + '\n'

            # get a dictionary for concept-ids
            concept2id = dict(
                zip(self.concepts, [i + 1 for i in range(len(self.concepts))]))
            for concept in self.concepts:
                out += '\n'
                indices = self.get_list(row=concept, flat=True)
                if self._modify_ref:
                    cogids = [self._modify_ref(self[i, ref]) for i in indices]
                else:
                    cogids = [self[i, ref] for i in indices]
                cogids = sorted(set(cogids))
                for cogid in cogids:
                    if cogid in self.msa[ref]:
                        for i, alm in enumerate(self.msa[ref][cogid]['alignment']):
                            taxon = self.msa[ref][cogid]['taxa'][i]
                            cid = concept2id[concept]
                            # add this line for alignments containing loans
                            real_cogid = self[self.msa[ref][cogid]['ID'][i], ref]
                            if not kw['confidence']:
                                alm_string = '\t'.join(alm)
                            else:
                                confs = [
                                    '{0}'.format(x)
                                    for x in self.msa[ref][cogid]['confidence'][i]]
                                chars = [
                                    x for x in self.msa[ref][cogid]['_charmat'][i]]
                                alm_string = '\t'.join(
                                    [a + '/' + b + '/' + c
                                     for a, b, c in zip(alm, confs, chars)]
                                )

                            out += util.tabjoin(
                                real_cogid, taxon, concept, cid, alm_string) + '\n'
                    else:
                        this_idx = [x for x in self.etd[ref][cogid] if x != 0][0][0]
                        taxon = self[this_idx, 'taxon']
                        seq = self[this_idx, 'ipa']
                        if not seq:
                            seq = ' '.join(self[this_idx, 'tokens'])
                        cid = concept2id[concept]
                        out += util.tabjoin(
                            cogid, taxon, concept, cid, ''.join(seq)) + '\n'

            util.write_text_file(filename + '.' + fileformat, out)

        if fileformat == 'msa':
            for key, value in sorted(self.msa[kw['ref']].items(), key=lambda x: x[0]):
                util.write_text_file(
                    os.path.join(
                        '{0}-msa'.format(value['dataset']),
                        '{0}-{1}.msa'.format(value['dataset'], key)),
                    msa2str(value, wordlist=kw['style'] in ['id', 'with_id']),
                    log=False)


def SCA(infile, **keywords):
    """
    Method returns alignment objects depending on input file or input data.

    Notes
    -----
    This method checks for the type of an alignment object and returns an
    alignment object of the respective type.
    """
    util.setdefaults(
        keywords,
        comment=rcParams['comment'],  # '#',
        diacritics=rcParams['diacritics'],  # None,
        vowels=rcParams['vowels'],  # None,
        tones=rcParams['tones'],  # None,
        combiners=rcParams['combiners'],  # '\u0361\u035c',
        breaks=rcParams['breaks'],  # '.-',
        stress=rcParams['stress'],  # "ˈˌ'",
        merge_vowels=rcParams['merge_vowels'],  # True
    )

    # check for datatype
    if isinstance(infile, dict):
        parent = MSA(infile, **keywords)
    else:
        cls_map = dict(msa=MSA, msq=MSA, psq=PSA, psa=PSA, csv=Alignments, tsv=Alignments)
        # lookup class by file extension:
        parent = cls_map[os.path.splitext(infile)[1][-3:]](infile, **keywords)

    return parent


def get_consensus(msa, gaps=False, taxa=False, classes=False, **keywords):
    """
    Calculate a consensus string of a given MSA.

    Parameters
    ----------
    msa : {c{list} ~lingpy.align.multiple.Multiple}
        Either an MSA object or an MSA matrix.
    gaps : c{bool} (default=False)
        If set to c{True}, return the gap positions in the consensus.
    taxa : {c{list} bool} (default=False)
        If *tree* is chosen as a parameter, specify the taxa in order of the aligned
        strings.
    classes : c{bool} (default=False)
        Specify whether sound classes shall be used to calculate the consensus.
    model : ~lingpy.data.model.Model
        A sound class model according to which the IPA strings shall be
        converted to sound-class strings.
    local : { c{bool}, "peaks", "gaps" }(default=False)
        Specify whether local pre-processing should be applied to the data. If
        set to c{peaks}, the average alignment score of each column is taken as
        reference to remove low-scoring columns from the alignment. If set to
        "gaps", the columns with the highest proportion of gaps will be
        excluded.

    Returns
    -------
    cons : c{str}
        A consensus string of the given MSA.
    """

    util.setdefaults(
        keywords,
        model=rcParams['sca'],
        gap_scale=1.0,
        mode='majority',
        gap_score=-10,
        weights=[1 for i in range(len(msa[0]))],
        local=False)

    # transform the matrix
    matrix = misc.transpose(getattr(msa, 'alm_matrix', msa))

    # check for local peaks
    if keywords['local']:
        if keywords['local'] == 'peaks':
            # calculate a local index
            peaks = []
            for line in matrix:
                sim = []
                for (i, charA), (j, charB) in util.combinations2(enumerate(line)):
                    if charA not in rcParams['gap_symbol'] and charB not in rcParams['gap_symbol']:
                        sim.append(keywords['model'](
                            tokens2class([charA], keywords['model'])[0],
                            tokens2class([charB], keywords['model'])[0]
                        ))
                    else:
                        sim.append(0.0)
                peaks.append(sum(sim) / len(sim))

            # get the average,min, and max of the peaks
            pmean = sum(peaks) / len(peaks)
            pmax = max(peaks)

            # exclude those lines from matrix whose average is smaller than pmean
            i = len(matrix) - 1
            for peak in peaks[::-1]:
                if peak <= pmax - pmean - pmean / 3:
                    del matrix[i]
                i -= 1
        elif keywords['local'] == 'gaps':
            # store the number of gaps in a simple array
            gap_array = []

            for line in matrix:
                gap_array += [line.count('-') / len(line)]

            # we now try to get the average number of lines
            average = sum(gap_array) / len(gap_array)

            # we discard all lines which are beyond the half of the average (stupid
            # solution, but for testing it hopefully suffices...)
            i = len(matrix) - 1
            for score in gap_array[::-1]:
                if score >= 1 - average + average / 4:
                    del matrix[i]
                i -= 1

    # check for classes
    if classes:
        # if classes are passed as array, we use this array as is
        if isinstance(classes, list):
            pass
        # if classes is a Model-object
        elif hasattr(msa, 'ipa2cls'):
            msa.ipa2cls(model=keywords['model'])
            classes = misc.transpose(msa.classes)
    # if no tree is passed, it is a simple majority-rule principle that outputs
    # the consensus string
    cons = []
    if not classes:
        for col in matrix:
            count = Counter(col)
            if rcParams['gap_symbol'] in count:
                count[rcParams['gap_symbol']] *= keywords['gap_scale']
            cons.append(count.most_common(1)[0][0])
    elif classes:
        for i, col in enumerate(classes):
            tmpA, tmpB = Counter(col), defaultdict(int)
            for j, c in enumerate(col):
                tmpB[matrix[i][j]] += 1
            # half the weight of gaps
            if rcParams['gap_symbol'] in tmpA:
                tmpA[rcParams['gap_symbol']] *= keywords['gap_scale']
            chars = tmpA.most_common()
            # if mode is set to 'maximize', calculate the score
            if keywords['mode'] == 'maximize':
                tmpC = {}
                for j, c in enumerate(col):
                    if c not in tmpC:
                        score = 0
                        for k, c2 in enumerate(col):
                            if '-' not in (c, c2):
                                score += keywords['model'](c, c2)
                            else:
                                if (c, c2) == (rcParams['gap_symbol'],
                                        rcParams['gap_symbol']):
                                    score += 0
                                else:
                                    score += keywords['gap_score'] * \
                                        keywords['weights'][i]

                        tmpC[c] = score / len(col)
                chars = [(c, n) for c, n in sorted(
                    tmpC.items(), key=lambda x:(x[1], tmpA[x[0]]), reverse=True)]
            # check for identical classes
            maxV = chars.pop(0)
            clss = [maxV[0]]
            while chars:
                newV = chars.pop(0)
                if newV[1] == maxV[1]:
                    clss += [newV[0]]
                else:
                    break
            tmp = Counter()
            for j, c in enumerate(col):
                if c in clss:
                    tmp.update([matrix[i][j]])
            chars = [c for c, n in tmp.most_common(2)]

            # apply check for gaps here, if there are more gaps than in the
            # full column, take the gaps, otherwise, take the next char
            if chars[0] == rcParams['gap_symbol']:
                if tmp[rcParams['gap_symbol']] > sum([tmp[x] for x in tmp if x != \
                    rcParams['gap_symbol']]):
                    cchar = rcParams['gap_symbol']
                else:
                    cchar = chars[1]
            else:
                cchar = chars[0]
            cons.append(cchar)

    return cons if gaps else [c for c in cons if c != rcParams['gap_symbol']]  
