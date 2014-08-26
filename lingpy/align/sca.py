# author   : Johann-Mattis List, Johannes Dellert
# email    : mattis.list@uni-marburg.de
# created  : 2013-03-07 20:07
# modified : 2014-07-22 13:46

"""
Basic module for pairwise and multiple sequence comparison.

The module consists of four classes which deal with pairwise and multiple
sequence comparison from the *sequence* and the *alignment* perspective. The
sequence perspective deals with unaligned sequences. The *alignment*
perspective deals with aligned sequences.

"""

__author__="Johann-Mattis List, Johannes Dellert"
__date__="2014-07-22"

import numpy as np
import re
import codecs
import os
import sys

from ..read.qlc import read_msa, normalize_alignment
from ..settings import rcParams
from ..basic.wordlist import Wordlist
from ..convert import html
from ..convert.tree import subGuideTree
from ..convert.strings import msa2str
from ..sequence.sound_classes import ipa2tokens, tokens2class, class2tokens, \
        prosodic_string, prosodic_weights
from .multiple import Multiple
from .pairwise import Pairwise
from ..algorithm import misc
from ._align import confidence

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

    def __init__(
            self, 
            infile, 
            **keywords
            ):

        # set the defaults
        defaults = {
                'comment'      : rcParams['comment'],
                "diacritics"   : rcParams['diacritics'],
                "vowels"       : rcParams['vowels'],
                "tones"        : rcParams['tones'],
                "combiners"    : rcParams['combiners'],
                "breaks"       : rcParams['breaks'],
                "stress"       : rcParams['stress'],
                "merge_vowels" : rcParams['merge_vowels'],
                "ids"          : False,
                "header"       : True,
                "normalize"    : True,
                }
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]
        
        # store comment-string
        self.comment = keywords['comment']

        # initialization checks first, whether we are dealing with msa-files or
        # with other, unaligned, sequence files and starts the
        # loading-procedures accordingly
        if type(infile) == dict:
            self._init_dict(infile,**keywords)
        else:
            if infile.endswith('.msa') or infile.endswith('.msq'):
                tmp = read_msa(infile,**keywords)
                self._init_dict(tmp,**keywords)

    def _init_dict(
            self,
            initdict,
            **keywords
            ):
        """
        Initialize by passing a dictionary with the relevant values.
        """
        
        for key in initdict:

            if key not in 'seqs':
                setattr(self,key,initdict[key])

        Multiple.__init__(
                self,
                initdict['seqs'],
                **keywords
                )

        if 'alignment' in initdict:
            self.alm_matrix = initdict['alignment']
        if 'local' in initdict:
            self.local = initdict['local']
        if 'swaps' in initdict:
            self.swaps = initdict['swaps']

    def ipa2cls(
            self,
            **keywords
            ):
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
        defaults = dict(
                model = rcParams['sca'],
                stress = rcParams['stress']
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        self.classes = []
        
        self.model = keywords['model']

        # redefine the sequences of the Multiple class
        class_strings = [tokens2class(
            seq.split(' '),
            self.model,
            stress=keywords['stress']) for seq in self.seqs]
        
        # define the scoring dictionaries according to the methods
        aligned_seqs = [alm for alm in self.alm_matrix]
        for i in range(len(aligned_seqs)):
            self.classes.append(
                    list(
                        ''.join(
                            class2tokens(
                                class_strings[i],
                                aligned_seqs[i]
                                )
                            ).replace('-','X')
                        )
                    )
    def output(
            self,
            fileformat = 'msa',
            filename =  None,
            sorted_seqs = False,
            unique_seqs = False,
            **keywords
            ):
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
        defaults = dict(
                wordlist = False,
                timestamp = False
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]
                
        if not filename:
            filename = self.infile

        # define the outfile and check, whether it already exists
        outfile = filename + '.' + fileformat

        # create a specific format string in order to receive taxa of equal
        # length
        mtax = max([len(t) for t in self.taxa])
        txf = '{0:.<'+str(mtax)+'}'
        
        if fileformat not in ['html', 'tex']:
            out = codecs.open(outfile,'w','utf-8')

            # start writing data to file
            out.write(self.dataset+'\n')

        if fileformat in ['msq','msa']:
            out.write(self.seq_id+'\n')

        if not sorted_seqs or fileformat == 'psa':
            for i,taxon in enumerate(self.taxa):
                if fileformat == 'msq':
                    out.write(txf.format(taxon)+'\t'+self.seqs[i]+'\n')
                elif fileformat == 'msa':
                    out.write(txf.format(taxon)+'\t')
                    out.write('\t'.join(self.alm_matrix[i])+'\n')
                elif fileformat == 'psa':
                    if not hasattr(self,'alignments'):
                        self.get_pairwise_alignments(new_calc=False)
                    else:
                        pass
                    for j,taxonB in enumerate(self.taxa):
                        if i < j:
                            try:
                                almA,almB,score = self.alignments[i,j]
                            except:
                                almB,almA,score = self.alignments[j,i]
                            
                            out.write('{0} ({1}, {2})\n'.format(
                                self.seq_id,
                                taxon,
                                taxonB))
                            out.write(txf.format(taxon)+'\t')
                            out.write('\t'.join(almA)+'\n')
                            out.write(txf.format(taxonB)+'\t')
                            out.write('\t'.join(almB)+'\n')
                            out.write('{0} {1:.2f}\n\n'.format(
                                self.comment,
                                score))
        elif sorted_seqs:
            if fileformat == 'msa':
                alms = ['\t'.join(alm) for alm in self.alm_matrix]
            else:
                alms = [seq for seq in self.seqs]
            
            if not unique_seqs:
                taxalms = zip(self.taxa,alms)
                taxalms = sorted(taxalms,key=lambda x:x[1])
            
            elif unique_seqs:
                uniqs = sorted([x[0] for x in self.uniseqs.values()])
                taxa = [self.taxa[x] for x in uniqs]
                alms = [alms[x] for x in uniqs]
                taxalms = zip(taxa,alms)
                taxalms = sorted(taxalms,key=lambda x:x[1])

            for taxon,alm in taxalms:
                out.write(txf.format(taxon)+'\t'+alm+'\n')


        if fileformat == 'msa':
            if hasattr(self,'local'):
                if self.local:
                    out.write(txf.format("LOCAL")+'\t')
                    tmp = ['.'] * len(self.alm_matrix[0])
                    for i in self.local:
                        tmp[i] = '*'
                    out.write('\t'.join(tmp)+'\n')

            if hasattr(self,'swaps'):
                if self.swaps:
                    out.write(txf.format('SWAPS')+'\t')
                    tmp = ['.'] * len(self.alm_matrix[0])
                    for i in self.swaps:
                        tmp[i[0]] = '+'
                        tmp[i[1]] = '-'
                        tmp[i[2]] = '+'
                    out.write('\t'.join(tmp)+'\n')
            if hasattr(self,'merge'):
                if len(set(self.merge.values())) < len(set(self.merge.keys())):
                    out.write(txf.format('MERGE')+'\t')
                    tmp = ['.'] * len(self.alm_matrix[0])
                    start = False
                    before = 1
                    for k in sorted(self.merge):
                        if self.merge[k] == before and not start:
                            tmp[k-1] = '<'
                            start = True
                        if self.merge[k] == before and start:
                            tmp[k] = '-'
                        if self.merge[k] != before and start:
                            start = False
                            tmp[k-1] = '>'
                            before += 1
                        elif self.merge[k] != before and not start:
                            before += 1
                    out.write('\t'.join(tmp)+'\n')
            if hasattr(self,'proto'):
                out.write(txf.format('PROTO')+'\t')
                out.write('\t'.join(self.proto)+'\n')
            if hasattr(self,'consensus'):
                out.write(txf.format("CONSE")+'\t')
                out.write('\t'.join(self.consensus)+'\n')



        if fileformat in ['html','tex']:
            self.output('msa', '.tmp', sorted_seqs, unique_seqs)
            if 'filename' not in keywords:
                keywords['input_file'] = os.path.split(self.infile)[1]
                keywords['filename'] = filename

            if fileformat == 'html':
                html.msa2html('.tmp.msa', **keywords)
            elif fileformat == 'tex':
                html.msa2tex('.tmp.msa', **keywords)
            try:
                os.remove('.tmp.msa')
            except:
                pass
        
        if fileformat not in ['html','tex'] and keywords['timestamp']:
            try:
                out.write('# Created using LingPy-2.2\n')
                if hasattr(self,'params'):
                    out.write('# Parameters: '+self.params+'\n')
                out.write('# Created: {0}\n'.format(rcParams['timestamp']))
            except:
                pass
            out.close()
        elif fileformat not in ['html', 'tex']:
            out.close()

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
    def __init__(
            self,
            infile,
            **keywords
            ):
        
        # set the defaults
        defaults = {
                'comment':rcParams['comment'],
                "diacritics" : rcParams['diacritics'],
                "vowels":rcParams['vowels'],
                "tones":rcParams['tones'],
                "combiners":rcParams['combiners'],
                "breaks":rcParams['breaks'],
                "stress":rcParams["stress"],
                "merge_vowels" : rcParams['merge_vowels']
                }

        # check for keywords
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        # add comment-char
        self.comment = keywords['comment']

        # check the ending of the infile
        if infile.endswith('.psa'):
            self._init_psa(infile,**keywords)
        else:
            self._init_seq(infile,**keywords)

    def _init_psa(
            self,
            infile,
            **keywords
            ):
        """
        Load a ``psa``-file.
        """
        # import the data from the input file
        self.infile = infile.split('/')[-1].replace('.psa','')

        data = []
        try:
            raw_data = codecs.open(infile+'.psa','r','utf-8')
        except:
            raw_data = codecs.open(infile,'r','utf-8')

        for line in raw_data:
            if not line.startswith(self.comment):
                data.append(line.strip())

        # set the first parameters
        self.dataset = data[0]
        
        # delete the first line of the data, since they are no longer needed
        del data[0]

        # append the other lines of the data, they consist of triplets,
        # separated by double line breaks
        self.taxa = []
        self.pairs = []
        self.seq_ids = []
        self.alignments = []

        i = 0
        while i <= len(data) - 3:
            try:
                self.seq_ids.append(data[i])
                
                datA = data[i+1].split('\t')
                datB = data[i+2].split('\t')
                
                taxonA = datA[0]
                taxonB = datB[0]
                almA = datA[1:]
                almB = datB[1:]
                
                self.taxa.append((taxonA,taxonB))
                self.pairs.append(
                        (
                            '.'.join([k for k in almA if k != '-']),
                            '.'.join([k for k in almB if k != '-'])
                            )
                        )
                self.alignments.append(
                        (
                            [str(a) for a in almA],
                            [str(b) for b in almB],
                            0)
                        )
                i += 4
            except:
                print("[!] Line {0} of the data is probablyb miscoded.".format(
                    i+1
                    ))
                i += 1

        Pairwise.__init__(
                self,
                self.pairs,
                **keywords
                )

    def _init_seq(
            self,
            infile,
            **keywords
            ):
        """
        Load a ``psq``-file.
        """
        # import the data from the input file
        self.infile = infile.split('/')[-1].replace('.psq','')

        data = []
        try:
            raw_data = codecs.open(infile+'.psq','r','utf-8')
        except:
            raw_data = codecs.open(infile,'r','utf-8')

        for line in raw_data:
            if not line.startswith(self.comment):
                data.append(line.strip())

        # set the first parameters
        self.dataset = data[0]
        
        # delete the first line of the data, since they are no longer needed
        del data[0]

        # append the other lines of the data, they consist of triplets,
        # separated by double line breaks
        self.taxa = []
        self.pairs = []
        self.seq_ids = []
        i = 0
        while i <= len(data) - 3:
            try:
                self.seq_ids.append(data[i])
                taxonA,seqA = data[i+1].split('\t')
                taxonB,seqB = data[i+2].split('\t')
                self.taxa.append((taxonA.strip('.'),taxonB.strip('.')))
                self.pairs.append((seqA,seqB))
                i += 4
            except:
                print("[!] Line "+str(i+1)+" of the data is probably miscoded.")

        self.pair_num = len(self.pairs)

        Pairwise.__init__(
                self,
                self.pairs,
                **keywords
                )

    def output(
            self,
            fileformat = 'psa',
            filename =  None,
            **keywords
            ):
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
        defaults = dict(
                gop = -2,
                model = rcParams['sca'],
                transform = rcParams['align_transform'],
                scores = False
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        if not filename:
            filename = self.infile

        # define the outfile and check, whether it already exists
        outfile = filename + '.' + fileformat
        # check whether outfile already exists
        try:
            tmp = codecs.open(outfile,'r','utf-8')
            tmp.close()
            outfile = filename + '_out.' + fileformat
        except:
            pass

        # open output file
        out = codecs.open(outfile,'w','utf-8')

        # if data is simple, just write simple data to file
        if fileformat == 'psq':
            out.write(self.dataset + '\n')
            for i,(a,b) in enumerate(self.pairs):
                out.write(self.seq_ids[i]+'\n')
                
                # determine longest taxon in order to create a format string
                # for taxa of equal length
                mtax = max([len(t) for t in self.taxa[i]])
                txf = '{0:.<'+str(mtax)+'}'

                out.write(txf.format(self.taxa[i][0])+'\t'+a+'\n')
                out.write(txf.format(self.taxa[i][1])+'\t'+b+'\n\n')
            out.close()

        # if data is psa-format
        elif fileformat == 'psa':
            out.write(self.dataset + '\n')
            for i,(a,b,c) in enumerate(self.alignments):
                out.write(self.seq_ids[i]+'\n')
                
                # determine longest taxon in order to create a format string
                # for taxa of equal length
                mtax = max([len(t) for t in self.taxa[i]])
                txf = '{0:.<'+str(mtax)+'}'
                
                out.write(txf.format(self.taxa[i][0])+'\t'+'\t'.join(a)+'\n')
                out.write(txf.format(self.taxa[i][1])+'\t'+'\t'.join(b)+'\n')
                
                if keywords['scores']:
                    # get partial alignment scores
                    scores = []
                    idxA,idxB = 0,0
                    proA = self.weights[i][0] 
                    proB = self.weights[i][1]

                    for x,y in zip(a,b):
                        if '-' not in (x,y):
                            try:
                                scores += [self.model(x,y)]
                            except:
                                self._set_model(model=keywords['model'])
                            idxA += 1
                            idxB += 1

                        else:
                            if x == '-':
                                scores += [keywords['gop'] * proB[idxB]]
                                idxB += 1
                            elif y == '-':
                                scores += [keywords['gop'] * proA[idxA]]
                                idxA += 1

                            
                    out.write(
                            txf.format(self.comment)+'\t'+'\t'.join(
                                ['{0:.2f}'.format(s) for s in scores]
                                )+'\t{0:.2f}\n'.format(sum(scores))
                            )
                out.write('{0} {1:.2f}'.format(self.comment,c)+'\n\n')
            out.close()

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
    loans : bool (default=True)
        Specify whether loans should be included in the cognate sets.

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
            row = 'concept',
            col = 'doculect',
            conf = '',
            ref = 'cogid', 
            loans = True,
            **keywords
            ):
        # initialize the wordlist
        Wordlist.__init__(self,infile,row,col,conf)
        
        # check for reference / cognates
        if 'cognates' in keywords:
            print(rcParams['W_deprecation'].format('cognates','ref'))
            ref = keywords['cognates']
        
        # change ref to rcParams
        if ref != rcParams['ref']:
            rcParams['ref'] = ref

        # check for cognate-id or alignment-id in header
        try:
            self.etd = {ref:self.get_etymdict(ref=ref,loans=loans)}
        # else raise error
        except: 
            raise ValueError(
                    "[i] Did not find a cognate ID in the input file."
                    )

        # store loan-status
        self._loans = loans

        # check for strings
        if 'tokens' in self.header:
            stridx = self.header['tokens']
        elif 'orthoparse' in self.header:
            stridx = self.header['orthoparse']
        elif 'ipa' in self.header:
            stridx = self.header['ipa']
        else:
            if 'strings' in keywords:
                try:
                    stridx = self.header[keywords['strings']]
                except:
                    print(
                            "[i] No valid source for strings could be found."
                            )
            else:
                print("[i] No valid source for strings could be found.")
                return 

        # create the alignments by assembling the ids of all sequences
        if 'msa' not in self._meta:
            self._meta['msa'] = {ref:{}}
            for key,value in self.etd[ref].items():
                tmp = [x for x in value if x != 0]
                seqids = []
                for t in tmp:
                    seqids += t
                if len(seqids) > 1:
                    
                    # set up the dictionary
                    d = {}
                    d['taxa'] = []
                    d['seqs'] = []
                    d['dataset'] = self.filename
                    d['ID'] = []
                    d['alignment'] = []
                    if 'concept' in self.header:
                        concept = self[seqids[0],'concept']
                        d['seq_id'] = '{0} ("{1}")'.format(key,concept)
                    else:
                        d['seq_id'] = '{0}'.format(key)
                    
                    # set up the data
                    for seq in seqids:
                        taxon = self[seq,'taxa']
                        string = self[seq][stridx]
                        d['ID'] += [seq]
                        d['taxa'] += [self[seq,'taxa']]
                        d['seqs'] += [self[seq][stridx]]
                        if 'alignment' in self.header:
                            d['alignment'] += [self[seq,'alignment']]
                        else:
                            if type(string) == str:
                                d['alignment'] += [string.split(' ')]
                            else:
                                d['alignment'] += [string]
                                
                        d['alignment'] = normalize_alignment(d['alignment'])
                        
                    self._meta['msa'][ref][key] = d

    def _msa2col(
            self,
            ref='cogid'
            ):
        """
        Add alignments to column (space-separated) in order to make it easy to
        parse them in the wordlist editor.
        """
        tmp = {}
        for key,msa in self.msa[ref].items():
            for i,idx in enumerate(msa['ID']):
                try:
                    tmp[idx] = ' '.join(msa['alignment'][i])
                except KeyError:
                    print("[!] There are no alignments in your data.  Aborting...")
                    return
        missing = [idx for idx in self if idx not in tmp]
        for m in missing:
            if 'tokens' in self.header:
                tmp[m] = self[m,'tokens']
            elif 'ipa' in self.header:
                tmp[m] = ipa2tokens(self[m,'ipa'])
            elif 'alignment' in self.header:
                tmp[m] = self[m,'alignment']
            else:
                raise ValueError("There are no phonetic sequences (TOKENS, ALIGNMENT, or IPA) in your data.")

        self.add_entries('alignment', tmp, lambda x: x)

    def align(
            self,
            **keywords
            ):
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
        output : bool (default=False)
            Set to c{True} in order to write all alignments to file.
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
                method           = 'progressive',
                iteration        = False,
                swap_check       = False,
                output           = False,
                model            = rcParams['sca'],
                mode             = rcParams['align_mode'],
                modes            = rcParams['align_modes'],
                gop              = rcParams['align_gop'],
                scale            = rcParams['align_scale'],
                factor           = rcParams['align_factor'],
                tree_calc        = rcParams['align_tree_calc'],
                gap_weight       = rcParams['gap_weight'],
                restricted_chars = rcParams['restricted_chars'],
                classes          = rcParams['classes'],
                sonar            = rcParams['sonar'],
                scoredict        = rcParams['scorer'],
                ref              = rcParams['ref'],
                plots            = False,
                filename         = self.filename,
                show             = False,
                style            = 'plain',
                defaults         = False,
                )

        kw.update(keywords)
        if kw['defaults']: return kw

        if str(kw['model']) == kw['model']:
            kw['model'] = rcParams[kw['model']]

        # create a params attribute
        params = '_'.join(
                [
                    kw['method'],
                    kw['model'].name,
                    str(kw['gop']),
                    '{0:.1f}'.format(kw['scale']),
                    '{0:.1f}'.format(kw['factor']),
                    kw['tree_calc'],
                    '{0:.1f}'.format(kw['gap_weight']),
                    kw['restricted_chars']
                    ]
                )

        # define a score-bar that shows how far the work is processed
        if not rcParams['verbose'] and rcParams['_sverb']:
            task_len = len(self.msa[kw['ref']])
            if task_len >= rcParams['_sverb_tbar_len']:
                task_char = rcParams['_sverb_tchar']
                task_step = task_len / rcParams['_sverb_tbar_len']
                task_range = [int(x+0.5) for x in np.arange(0, task_len, task_step)] 
            else:
                task_range = list(range(task_len))
                task_char = rcParams['_sverb_tchar'] * int(rcParams['_sverb_tbar_len'] / task_len+0.5)
            task_string = ' ALIGNMENTS '.center(
                    rcParams['_sverb_tbar_len'],
                    rcParams['_sverb_fchar']
                    )
            task_string = '|' + task_string + '|'
            sys.stdout.write(task_string+'\r|')
            task_count = 0
            control_char = 0

        for key,value in sorted(
                self.msa[kw['ref']].items(),
                key=lambda x:x[0]
                ):
            if rcParams['verbose']: print("[i] Analyzing cognate set number {0}.".format(key))
            elif rcParams['_sverb']:
                if task_count in task_range and control_char < rcParams['_sverb_tbar_len']:
                    sys.stdout.write(task_char)
                    sys.stdout.flush()
                    control_char += len(task_char)
                task_count += 1
            
            # check for scorer keyword
            if not kw['scoredict']:
                m = SCA(
                        value,
                        **kw
                        )
            else:
                # get the tokens 
                numbers = [self[idx,'numbers'] for idx in value['ID']]
                if kw['sonar']:
                    sonars = [self[idx,'sonars'] for idx in value['ID']]
                else:
                    sonars = False
                value['seqs'] = numbers
                m = SCA(
                        value,
                        **kw
                        )

                kw['sonars'] = sonars
                kw['classes'] = False
            
            if kw['method'] == 'progressive':
                m.prog_align(**kw)
                        
            elif kw['method'] == 'library':
                m.lib_align(**kw)

            if kw['iteration']:
                m.iterate_clusters(0.5)
                m.iterate_orphans()
                m.iterate_similar_gap_sites()

            if kw['swap_check']:
                m.swap_check()

            # convert back to external format, if scoredict is set
            if kw['scoredict']:
                for i,alm in enumerate(m.alm_matrix):
                    tk = self[m.ID[i],'tokens']
                    new_tk = class2tokens(tk,alm)
                    m.alm_matrix[i] = new_tk

            self._meta['msa'][kw['ref']][key]['alignment'] = m.alm_matrix
            self._meta['msa'][kw['ref']][key]['_sonority_consensus'] = m._sonority_consensus
            self._meta['msa'][kw['ref']][key]['stamp'] = rcParams['align_stamp'].format(
                    m.dataset,
                    m.seq_id,
                    rcParams['timestamp'],
                    params
                    )


            if kw['output']:
                if kw['style'] in ['plain', 'msa']:
                    try:                        
                        m.output(
                                'msa',
                                filename='{0}-msa/{1}-{2}'.format(
                                    self.filename,
                                    m.dataset,
                                    key
                                    )
                                )
                    except:
                        os.mkdir('{0}-msa'.format(self.filename))
                        m.output(
                                'msa',
                                filename='{0}-msa/{1}-{2}'.format(
                                    self.filename,
                                    m.dataset,
                                    key
                                    )
                                )
                elif kw['style'] in ['with_id', 'id']:
                    msa_string = msa2str(
                            self._meta['msa'][kw['ref']][key],
                            wordlist=True
                            )
                    try:
                        f = codecs.open(
                                '{0}-msa/{1}-{2}'.format(
                                    self.filename,
                                    m.dataset,
                                    key
                                    ),
                                'w',
                                'utf-8'
                                )
                        f.write(msa_string)
                        f.close()
                    except:
                        os.mkdir('{0}-msa'.format(self.filename))
                        f = codecs.open(
                                '{0}-msa/{1}-{2}'.format(
                                    self.filename,
                                    m.dataset,
                                    key
                                    ),
                                'w',
                                'utf-8'
                                )
                        f.write(msa_string)
                        f.close()
        
        self._msa2col(kw['ref'])
        
        if not rcParams['verbose'] and rcParams['_sverb']: 
            if control_char < rcParams['_sverb_tbar_len']:
                sys.stdout.write(
                        (rcParams['_sverb_tbar_len'] - control_char) * rcParams['_sverb_tchar']
                            )
            sys.stdout.write('|\r'+rcParams['_sverb_tbar_len'] * ' '+'     \r')
            sys.stdout.flush()


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
        if rcParams['verbose']:
            print("[i] Successfully calculated confidence values for alignments.")
        
        return corrs

    def __len__(self):
        return len(self.msa)

    def _plot(
            self,
            fileformat = 'html',
            **keywords
            ):
        """
        Make an HTML plot of the aligned data.
        """
        defaults = dict(
                title = 'LexStat - Automatic Cognate Judgments',
                shorttitle = "LexStat",
                dataset = self.filename,
                show = False,
                filename = self.filename,
                ref = rcParams['ref'],
                confidence = False
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        self.output('alm',ref=keywords['ref'],filename='.tmp',
                confidence=keywords['confidence'])
        html.alm2html(
                '.tmp.alm',
                **keywords
                )
        try:
            os.remove('.tmp.alm')
        except:
            pass

    def get_consensus(
            self,
            tree = False,
            gaps = False,
            classes = False,
            consensus = 'consensus',
            counterpart = 'ipa',
            weights = [],
            **keywords
            ):
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
        # determine defaults
        defaults = dict(
                model     = rcParams['sca'],
                gap_scale = 1.0,
                ref       = rcParams['ref'],
                )
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        # check for deprecated "cognates"
        if 'cognates' in keywords:
            print(rcParams['W_deprecation'].format('cognates','ref'))
            ref = keywords['cognates']

        # switch ref
        if keywords['ref'] != rcParams['ref']:
            rcParams['ref'] = keywords['ref']

        # reassing ref for convenience
        ref = keywords['ref']

        # check for existing alignments
        test = list(self.msa[ref].keys())[0]
        if 'alignment' not in self.msa[ref][test]:
            print("[!] No alignments could be found, You should carry out"
                    " an alignment analysis first!")
            return
        
        # define a score-bar that shows how far the work is processed
        if not rcParams['verbose'] and rcParams['_sverb']:
            task_len = len(self.msa[ref])
            if task_len >= rcParams['_sverb_tbar_len']:
                task_char = rcParams['_sverb_tchar']
                task_step = task_len / rcParams['_sverb_tbar_len']
                task_range = [int(x+0.5) for x in np.arange(0, task_len, task_step)] 
            else:
                task_range = list(range(task_len))
                task_char = rcParams['_sverb_tchar'] * int(rcParams['_sverb_tbar_len'] / task_len+0.5)
            task_string = ' CONSENSUS '.center(
                    rcParams['_sverb_tbar_len'],
                    rcParams['_sverb_fchar']
                    )
            task_string = '|' + task_string + '|'
            sys.stdout.write(task_string+'\r|')
            task_count = 0
            control_char = 0


        # go on with the analysis
        cons_dict = {}
        for cog in self.etd[ref]:

            
            if cog in self.msa[ref]:
                if rcParams['verbose']: print("[i] Analyzing cognate set number '{0}'...".format(cog))
                elif rcParams['_sverb']:
                    if task_count in task_range and control_char < rcParams['_sverb_tbar_len']:
                        sys.stdout.write(task_char)
                        sys.stdout.flush()
                        control_char += len(task_char)
                    task_count += 1

                
                # temporary solution for sound-class integration
                if classes == True:
                    classes = []
                    if weights:
                        keywords['weights'] = prosodic_weights(
                                prosodic_string(
                                    self.msa[ref][cog]['_sonority_consensus']
                                    )
                                )
                    else:
                        keywords['weights'] = [1.0 for i in range(len(self.msa[ref][cog]['alignment']))]

                    for alm in self.msa[ref][cog]['alignment']:
                        cls = [c for c in tokens2class(
                                alm,
                                keywords['model']
                                ) if c != '0']
                        cls = class2tokens(cls,alm)
                        classes += [cls]

                    cons = get_consensus(
                            self.msa[ref][cog]['alignment'],
                            classes = misc.transpose(classes),
                            tree = tree,
                            gaps = gaps,
                            taxa = [str(taxon.replace("(","").replace(")","")) for taxon in self.msa[ref][cog]['taxa']],
                            #taxa = self.msa[ref][cog]['taxa'],
                            **keywords
                            )
                    classes = True
                else:
                    cons = get_consensus(
                            self.msa[ref][cog]['alignment'],
                            classes = classes,
                            tree = tree,
                            gaps = gaps,
                            taxa = [str(taxon.replace("(","").replace(")","")) for taxon in self.msa[ref][cog]['taxa']],
                            #taxa = self.msa[ref][cog]['taxa'],
                            **keywords
                            )
                self.msa[ref][cog]["consensus"] = cons

            # if there's no msa for a given cognate set, this set is a
            # singleton
            else:
                cons = self[[k[0] for k in self.etd[ref][cog] if k != 0][0],counterpart]
            
            # add consensus to dictionary
            cons_dict[cog] = cons        

        if not rcParams['verbose'] and rcParams['_sverb']: 
            if control_char < rcParams['_sverb_tbar_len']:
                sys.stdout.write(
                        (rcParams['_sverb_tbar_len'] - control_char) * rcParams['_sverb_tchar']
                            )
            sys.stdout.write('|\r'+rcParams['_sverb_tbar_len'] * ' '+'     \r')
            sys.stdout.flush()
        
        # add the entries
        self.add_entries(
                consensus,
                ref,
                lambda x:cons_dict[x] 
                )

    def output(
            self,
            fileformat,
            **keywords
            ):
        """
        Write wordlist to file.

        Parameters
        ----------
        fileformat : {"qlc", "tre","nwk","dst", "taxa", "starling", "paps.nex", "paps.csv" "html"}
            The format that is written to file. This corresponds to the file
            extension, thus 'csv' creates a file in csv-format, 'dst' creates
            a file in Phylip-distance format, etc. 
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
        cognates : str
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
        
        """
        kw = dict(
                ref = rcParams['ref'],
                filename = rcParams['filename'],
                style = "id",
                defaults = False,
                confidence = False,
                )
        kw.update(keywords)
        if kw['defaults']: return kw

        # check for html fileformat
        if fileformat == 'html':
            self._plot(**keywords)

        # define two vars for convenience
        ref = kw['ref']
        filename = kw['filename']

        if 'cognates' in kw:
            print(rcParams['W_deprecation'].format('cognates','ref'))
            ref = kw['cognates']

        if ref != rcParams['ref']:
            rcParams['ref'] = ref

        if fileformat not in ['alm', 'msa']:
            return self._output(fileformat,**kw)
        
        if fileformat == 'alm':

            # define the string to which the stuff is written
            out = self.filename+'\n'
            
            # get a dictionary for concept-ids
            concept2id = dict(
                    zip(
                        self.concepts,
                        [i+1 for i in range(len(self.concepts))]
                            )
                        )
            idx = 1
            for concept in self.concepts:
                out += '\n'
                indices = self.get_list(
                        row=concept,
                        flat=True
                        )
                if self._loans:
                    cogids = [abs(self[i,ref]) for i in indices]
                else:
                    cogids = [self[i,ref] for i in indices]
                cogids = sorted(set(cogids))
                for cogid in cogids:
                    if cogid in self.msa[ref]:
                        for i,alm in enumerate(self.msa[ref][cogid]['alignment']):
                            taxon = self.msa[ref][cogid]['taxa'][i]
                            seq = self.msa[ref][cogid]['seqs'][i]
                            cid = concept2id[concept]
                            # add this line for alignments containing loans
                            real_cogid = self[self.msa[ref][cogid]['ID'][i],ref]
                            if not kw['confidence']:
                                alm_string = '\t'.join(alm)
                            else:
                                confs = ['{0}'.format(x) for x in
                                        self.msa[ref][cogid]['confidence'][i]]
                                chars = [x for x in
                                        self.msa[ref][cogid]['_charmat'][i]]
                                alm_string = '\t'.join(
                                        [a+'/'+b+'/'+c for a,b,c in zip(
                                            alm,
                                            confs,
                                            chars
                                            )
                                            ]
                                        )

                            out += '\t'.join(
                                [
                                    str(real_cogid),
                                    taxon,
                                    concept,
                                    str(cid),
                                    alm_string, #'\t'.join(alm)
                                    ]
                                )+'\n'


                    else:
                        this_idx = [x for x in self.etd[ref][cogid] if x != 0][0][0]
                        taxon = self[this_idx,'taxon']
                        seq = self[this_idx,'ipa']
                        if not seq:
                            seq = ' '.join(self[this_idx,'tokens'])
                        cid = concept2id[concept]
                        out += '\t'.join(
                                [
                                    str(cogid),
                                    taxon,
                                    concept,
                                    str(cid),
                                    ''.join(seq),
                                    ]
                                )+'\n'

            f = codecs.open(filename + '.' + fileformat,'w','utf-8')
            f.write(out)
            f.close()

        if fileformat == 'msa':
            if kw['style'] in ['id', 'with_id']:
                wordlist = True
            else:
                wordlist = False

            for key,value in sorted(
                    self.msa[kw['ref']].items(),
                    key=lambda x:x[0]
                    ):
                
                msa_string = msa2str(
                        value,
                        wordlist=wordlist
                        )
                try:
                    f = codecs.open(
                            '{0}-msa/{1}-{2}.msa'.format(
                                self.filename,
                                value['dataset'],
                                key
                                ),
                            'w',
                            'utf-8'
                            )
                    f.write(msa_string)
                    f.close()

                except:
                    os.mkdir('{0}-msa'.format(self.filename))
                    f = codecs.open(
                            '{0}-msa/{1}-{2}.msa'.format(
                                self.filename,
                                value['dataset'],
                                key
                                ),
                            'w',
                            'utf-8'
                            )
                    f.write(msa_string)
                    f.close()

def SCA(
        infile,
        **keywords
        ):
    """
    Method returns alignment objects depending on input file or input data.

    Notes
    -----
    This method checks for the type of an alignment object and returns an
    alignment object of the respective type.
    """

    # set the defaults
    defaults = {
            'comment'      : rcParams['comment'], #'#',
            "diacritics"   : rcParams['diacritics'], #None,
            "vowels"       : rcParams['vowels'], #None,
            "tones"        : rcParams['tones'], #None,
            "combiners"    : rcParams['combiners'], #'\u0361\u035c',
            "breaks"       : rcParams['breaks'], #'.-',
            "stress"       : rcParams['stress'], #"ˈˌ'",
            "merge_vowels" : rcParams['merge_vowels'], #True
            }
    # check for keywords
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]
   
    # check for datatype
    if type(infile) == dict:
        parent = MSA(infile,**keywords)
    else:
        if infile[-4:] in ['.msa','.msq']:
            parent = MSA(infile,**keywords)
        elif infile[-4:] in ['.psq','psa']:
            parent = PSA(infile,**keywords)
        elif infile[-4:] in ['.csv','.tsv']:
            parent = Alignments(infile,**keywords)

    return parent

def get_consensus(
        msa,
        tree = False,
        gaps = False,
        taxa = False,
        classes = False,
        **keywords
        ):
    """
    Calculate a consensus string of a given MSA.

    Parameters
    ----------
    msa : {c{list} ~lingpy.align.multiple.Multiple}
        Either an MSA object or an MSA matrix.
    tree : {c{str} ~lingpy.thirdparty.cogent.PhyloNode}
        A tree object or a Newick string along which the consensus shall be
        calculated.
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
    # set defaults
    def rep_weights(char1, char2, prosody):
                    if char1 == char2:
                        return 0
                    else:
                        if char1 == '-':
                            if char2 in ['a','e','i','o','u','E','I',
                                       '3', 'ɛ', 'æ', 'ɜ', 'ɐ', 'ʌ', 
                                       'ᴇ', 'ə', 'ɘ', 'ɤ', 'ᴀ', 'ã', 
                                       'ɑ', 'ɪ', 'ɨ', 'ɿ', 'ʅ', 'ɯ',
                                       'œ', 'ɞ', 'ɔ', 'ø', 'ɵ', 'õ', 
                                       'ɶ', 'ɷ', 'ʏ', 'ʉ', 'ᴜ', 'ʊ']:
                                return 0.1
                            else:
                                return 0.8
                        elif char2 == "-":
                            if char1 in ['a','e','i','o','u','E','I',
                                       '3', 'ɛ', 'æ', 'ɜ', 'ɐ', 'ʌ', 
                                       'ᴇ', 'ə', 'ɘ', 'ɤ', 'ᴀ', 'ã', 
                                       'ɑ', 'ɪ', 'ɨ', 'ɿ', 'ʅ', 'ɯ',
                                       'œ', 'ɞ', 'ɔ', 'ø', 'ɵ', 'õ', 
                                       'ɶ', 'ɷ', 'ʏ', 'ʉ', 'ᴜ', 'ʊ']:
                                return 0.3
                            else:
                                return 0.1
                        else:
                            return 0.1
    defaults = dict(
            model = rcParams['sca'],
            gap_scale = 1.0,
            mode = 'majority',
            gap_score = -10,
            weights = [1 for i in range(len(msa[0]))],
            rep_weights = rep_weights,
            local = False
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]
    
    # stores the consensus string
    cons = []

    # transform the matrix
    if hasattr(msa,'alm_matrix'):
        matrix = misc.transpose(msa.alm_matrix)
    else:
        matrix = misc.transpose(msa)

    # check for local peaks
    if keywords['local']:
        
        if keywords['local'] == 'peaks':
            # calculate a local index
            peaks = []
            for line in matrix:
                sim = []
                for i,charA in enumerate(line):
                    for j,charB in enumerate(line):
                        if i < j:
                            if charA not in '-' and charB not in '-':
                                sim += [keywords['model'](
                                        tokens2class([charA],keywords['model'])[0],
                                        tokens2class([charB],keywords['model'])[0]
                                        )]
                                a = tokens2class([charA],keywords['model'])[0]
                                b = tokens2class([charB],keywords['model'])[0]

                            else:
                                sim += [0.0]
                peaks += [sum(sim) / len(sim)]

            # get the average,min, and max of the peaks
            pmean = sum(peaks) / len(peaks)
            pmax = max(peaks)
            pmin = min(peaks)

            # exclude those lines from matrix whose average is smaller than pmean
            i = len(matrix)-1
            for peak in peaks[::-1]:
                if peak <= pmax - pmean - pmean/3:
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
                if score >= 1 - average + average/4:
                    del matrix[i]
                i -= 1


    # check for classes
    if classes:
        # if classes are passed as array, we use this array as is
        if type(classes) == list:
            pass
        # if classes is a Model-object
        elif hasattr(msa,'ipa2cls'):
            msa.ipa2cls(model=keywords['model'])
            classes = misc.transpose(msa.classes)
    
    # if no tree is passed, it is a simple majority-rule principle that outputs
    # the consensus string
    if not tree:
        if not classes:
            for col in matrix:
                tmp = {}

                # count chars in columns
                for j,c in enumerate(col):
                    if c in tmp:
                        tmp[c] += 1
                    else:
                        tmp[c] = 1

                # half the weight of gaps
                if '-' in tmp:
                    tmp['-'] = tmp['-'] * keywords['gap_scale']

                # get maximum
                chars = [c for c,n in sorted(
                    tmp.items(),
                    key=lambda x:(x[1],len(x[0])),
                    reverse=True
                    )]
                
                # append highest-scoring char
                cons += [chars[0]]
        elif classes:
            for i,col in enumerate(classes):
                tmpA = {}
                tmpB = {}

                # count chars in columns
                for j,c in enumerate(col):
                    if c in tmpA:
                        tmpA[c] += 1
                    else:
                        tmpA[c] = 1

                    if matrix[i][j] in tmpB:
                        tmpB[matrix[i][j]] += 1
                    else:
                        tmpB[matrix[i][j]] = 1

                # half the weight of gaps
                if '-' in tmpA:
                    tmpA['-'] = tmpA['-'] * keywords['gap_scale']
                    #print(tmp['-'],keywords['gap_scale'])               

                # get max
                chars = [(c,n) for c,n in sorted(
                    tmpA.items(),
                    key=lambda x:x[1],
                    reverse=True
                    )]

                # if mode is set to 'maximize', calculate the score 
                if keywords['mode'] == 'maximize':
                    tmpC = {}
                    for j,c in enumerate(col):
                        if c in tmpC:
                            pass
                        else:
                            score = 0
                            for k,c2 in enumerate(col):
                                if '-' not in (c,c2):
                                    score += keywords['model'](c,c2)
                                else:
                                    if (c,c2) == ('-','-'):
                                        score += 0
                                    else:
                                        score += keywords['gap_score'] * keywords['weights'][i]

                            score = score / len(col)
                            tmpC[c] = score

                    chars = [(c,n) for c,n in sorted(
                        tmpC.items(),
                        key=lambda x:(x[1],tmpA[x[0]]),
                        reverse=True
                        )]

                # check for identical classes
                maxV = chars.pop(0)

                clss = [maxV[0]]
                while chars:
                    newV = chars.pop(0)
                    if newV[1] == maxV[1]:
                        clss += [newV[0]]
                    else:
                        break
                
                tmp = {}
                for j,c in enumerate(col):
                    if c in clss:
                        if matrix[i][j] in tmp:
                            tmp[matrix[i][j]] += 1
                        else:
                            tmp[matrix[i][j]] = 1
 
                # get max
                chars = [c for c,n in sorted(
                    tmp.items(),
                    key=lambda x:(x[1],len(x[0])),
                    reverse=True
                    )]

                # apply check for gaps here, if there are more gaps than in the
                # full column, take the gaps, otherwise, take the next char
                if chars[0] == '-':
                    if tmp['-'] > sum([tmp[x] for x in tmp if x != '-']):
                        cchar = '-'
                    else:
                        cchar = chars[1]
                else:
                    cchar = chars[0]
                    
                cons += [cchar]
                #cons += [chars[0]]

    # otherwise, we use a bottom-up parsimony approach to determine the best
    # match
    elif tree:
        #hack: this operates on the non-transposed alm_matrix
        matrix = misc.transpose(matrix)
        #if the taxa are defined, extract the sub-guidetree accordingly
        if taxa:
            all_taxa = [leaf.Name for leaf in tree.tips()]
            taxon_to_id = dict({(str(all_taxa[i]),i) for i in range(0,len(all_taxa))})
            #print(taxon_to_id)
            tree = tree.deepcopy()
            for leaf in tree.tips():
                leaf.Name = str(taxon_to_id[leaf.Name])
            #print(tree)
            #print(taxa)
            newIndices = [taxon_to_id[str(taxon)] for taxon in taxa]
            #print(newIndices)
            tree = subGuideTree(tree,newIndices)
            #print(tree)
            #indexMap = dict({(newIndices[taxa[i]],i) for i in range(len(taxa))})
            #print(indexMap)
            #for leaf in tree.tips():
            #    leaf.Name = str(indexMap[int(leaf.Name)])
        #otherwise, the leaves of the guide trees are expected to have integer IDs
        #print("\nWrite partial alignments and sizes into the tree:")
        for node in tree.postorder():
            if node.isTip():
                #print(node.Name)
                node.alignment = [matrix[int(node.Name)]]
                node.size = 1
            else:
                node.alignment = node.Children[0].alignment + node.Children[1].alignment
                node.size = node.Children[0].size + node.Children[1].size
        #for node in tree.postorder():
        #    print str(node.Name) + ": (" + str(node.size) + ") " + str(node.alignment)
        #print("\nCompute phoneme distribution at each position of the alignment:")
        for node in tree.postorder():
            node.distribution = []
        for i in range (0,len(matrix[0])):
            for node in tree.postorder():
                if node.isTip():
                    node.distribution.append({matrix[int(node.Name)][i] : 1.0})
                else:
                    node.distribution.append({})
                    child1 = node.Children[0]
                    child2 = node.Children[1]
                    for phoneme in set(child1.distribution[i].keys()) | set(child2.distribution[i].keys()):
                        value = 0.0
                        if phoneme in child1.distribution[i].keys():
                            value += child1.size * child1.distribution[i][phoneme]
                        if phoneme in child2.distribution[i].keys():
                            value += child2.size * child2.distribution[i][phoneme]
                        value /= node.size
                        node.distribution[i][phoneme] = value
        #for node in tree.postorder():
        #    print str(node.Name) + ": " + str(node.distribution)
        recon_alg = "sankoff_parsimony"
        #print("\nReconstruct word forms at inner nodes by simplistic criteria:")
        for node in tree.postorder():
            node.reconstructed = []
        if recon_alg == "modified_consensus":
            for i in range (0,len(matrix[0])):
                for node in tree.postorder():
                    dist = node.distribution[i]
                    maxValue = max(dist.values())
                    maxKeys = [key for key in dist.keys() if dist[key]==maxValue]
                    if len(maxKeys) == 1 or node.isRoot():
                        node.reconstructed.append(maxKeys[0])
                    else:         
                        parentDist = node.Parent.distribution[i]
                        maxKey = max(maxKeys, key=(lambda key: parentDist[key]))
                        node.reconstructed.append(maxKey)
        elif recon_alg == "binary_decision":
            for i in range (0,len(matrix[0])):
                for node in tree.postorder():
                    if node.isTip():
                        #just retrieve the original at the leaves
                        dist = node.distribution[i]
                        maxValue = max(dist.values())
                        maxKeys = [key for key in dist.keys() if dist[key]==maxValue]
                        node.reconstructed.append(maxKeys[0])
                    else:
                        leftVariant = node.Children[0].reconstructed[i]
                        rightVariant = node.Children[1].reconstructed[i]
                        #if one of both is '-', take the other one (preference for segment loss)
                        #BUT: epenthesis is allowed
                        if leftVariant == '-' and rightVariant not in ['a','e','i','o','u','E','I',
                                                                       '3', 'ɛ', 'æ', 'ɜ', 'ɐ', 'ʌ', 
                                                                       'ᴇ', 'ə', 'ɘ', 'ɤ', 'ᴀ', 'ã', 
                                                                       'ɑ', 'ɪ', 'ɨ', 'ɿ', 'ʅ', 'ɯ',
                                                                       'œ', 'ɞ', 'ɔ', 'ø', 'ɵ', 'õ', 
                                                                       'ɶ', 'ɷ', 'ʏ', 'ʉ', 'ᴜ', 'ʊ']:
                            node.reconstructed.append(rightVariant)
                        elif rightVariant == '-' and leftVariant not in ['a','e','i','o','u','E','I',
                                                                       '3', 'ɛ', 'æ', 'ɜ', 'ɐ', 'ʌ', 
                                                                       'ᴇ', 'ə', 'ɘ', 'ɤ', 'ᴀ', 'ã', 
                                                                       'ɑ', 'ɪ', 'ɨ', 'ɿ', 'ʅ', 'ɯ',
                                                                       'œ', 'ɞ', 'ɔ', 'ø', 'ɵ', 'õ', 
                                                                       'ɶ', 'ɷ', 'ʏ', 'ʉ', 'ᴜ', 'ʊ']:
                            node.reconstructed.append(leftVariant)
                        else:
                            #let the distribution decide otherwise
                            dist = node.distribution[i]
                            maxValue = max(dist.values())
                            maxKeys = [key for key in dist.keys() if dist[key]==maxValue]
                            if len(maxKeys) == 1 or node.isRoot():
                                node.reconstructed.append(maxKeys[0])
                            else:         
                                parentDist = node.Parent.distribution[i]
                                maxKey = max(maxKeys, key=(lambda key: parentDist[key]))
                                node.reconstructed.append(maxKey)
        elif recon_alg == "sankoff_parsimony":
            systematizationBonus = True
            distributionBonus = True
            for node in tree.postorder():
                node.sankoffTable = []
                node.sankoffPointers = []
            for i in range (0,len(matrix[0])):
                def sankoff_value(sankoffTable, char):
                    if char in sankoffTable.keys():
                        return sankoffTable[char]
                    else:
                        return 65535; #approximation to Integer.MAX_INT
                mtx = keywords["rep_weights"]
                #apply the Sankoff algorithm, store backpointers
                for node in tree.postorder():
                    node.sankoffTable.append(dict())
                    dist = node.distribution[i]
                    if node.isTip():
                        #just retrieve the original at the leaves (needs to be reconstructed from the differences)
                        maxValue = max(dist.values())
                        maxKeys = [key for key in dist.keys() if dist[key]==maxValue]
                        node.sankoffTable[i][maxKeys[0]] = 0
                        if systematizationBonus:
                            node.orig.best_replacements = {}
                    else:
                        node.sankoffPointers.append(dict())
                        sankoff1 = node.Children[0].sankoffTable[i]
                        sankoff2 = node.Children[1].sankoffTable[i]
                        recon1 = node.Children[0].reconstructed
                        recon2 = node.Children[1].reconstructed
                        for char in set(sankoff1.keys()) | set(sankoff2.keys()):
                            minSankoffValue = 65536
                            backPointers = ['-','-']
                            for char1 in sankoff1.keys():
                                for char2 in sankoff2.keys():
                                    bonusFactor1 = 1.0
                                    bonusFactor2 = 1.0
                                    if systematizationBonus:
                                        if hasattr(node.Children[0].orig, "best_replacements") and char + char1 in node.Children[0].orig.best_replacements.keys():
                                            bonusFactor1 /= 1 + node.Children[0].orig.best_replacements[char+char1] 
                                            #print ("bonusFactor1: " + str(bonusFactor1))           
                                        if hasattr(node.Children[1].orig, "best_replacements") and char + char2 in node.Children[1].orig.best_replacements.keys():
                                            bonusFactor2 /= 1 + node.Children[1].orig.best_replacements[char+char2]
                                            #print ("bonusFactor2: " + str(bonusFactor2))   
                                    if distributionBonus:
                                        bonusFactor1 /= node.distribution[i][char]
                                        bonusFactor2 /= node.distribution[i][char]
                                    char1no5 = char1
                                    char2no5 = char2
                                    if char1no5 == "5": char1no5 = "N"
                                    if char2no5 == "5": char2no5 = "N"
                                    proso1 = prosodic_string(recon1 + [char1no5])[-1]
                                    proso2 = prosodic_string(recon2 + [char2no5])[-1]
                                    sankoffValue = mtx(char,char1,proso1) * bonusFactor1 + sankoff1[char1]  + mtx(char,char2,proso2) * bonusFactor2 + sankoff2[char2]
                                    if (sankoffValue < minSankoffValue):
                                        minSankoffValue = sankoffValue
                                        backPointers[0] = char1
                                        backPointers[1] = char2
                            node.sankoffTable[i][char] = minSankoffValue
                            node.sankoffPointers[i][char] = backPointers
                        if systematizationBonus:           
                            if not hasattr(node.Children[0].orig, "best_replacements"):
                               node.Children[0].orig.best_replacements = {}
                            if not hasattr(node.Children[1].orig, "best_replacements"):
                               node.Children[1].orig.best_replacements = {}
                            minValue = min(node.sankoffTable[i].values())
                            minKeys = [key for key in node.sankoffTable[i].keys() if node.sankoffTable[i][key]==minValue]
                            for char in minKeys:
                                minValue1 = min(node.Children[0].sankoffTable[i].values())
                                for char1 in [key for key in node.Children[0].sankoffTable[i].keys() if node.Children[0].sankoffTable[i][key]==minValue1]:
                                    if char + char1 not in node.Children[0].orig.best_replacements.keys():
                                        node.Children[0].orig.best_replacements[char + char1] = 1
                                    else:
                                        node.Children[0].orig.best_replacements[char + char1] += 1
                                    #print(node.Children[0].orig.best_replacements)
                                minValue2 = min(node.Children[1].sankoffTable[i].values())
                                for char2 in [key for key in node.Children[1].sankoffTable[i].keys() if node.Children[1].sankoffTable[i][key]==minValue2]:
                                    if char + char2 not in node.Children[1].orig.best_replacements.keys():
                                        node.Children[1].orig.best_replacements[char + char2] = 1
                                    else:
                                        node.Children[1].orig.best_replacements[char + char2] += 1
                                    #print(node.Children[1].orig.best_replacements)
                #read out the backpointers for an optimal reconstruction
                #print(str(tree.sankoffTable[i]))
                minValue = min(tree.sankoffTable[i].values())
                minKeys = [key for key in tree.sankoffTable[i].keys() if tree.sankoffTable[i][key]==minValue]
                reconChar = minKeys[0]
                def reconstruct_by_sankoff(node, char):
                    if char == "5": node.reconstructed.append("N")
                    else: node.reconstructed.append(char)
                    if not node.isTip():
                        reconstruct_by_sankoff(node.Children[0], node.sankoffPointers[i][char][0])
                        reconstruct_by_sankoff(node.Children[1], node.sankoffPointers[i][char][1])
                reconstruct_by_sankoff(tree, reconChar)
        else:
            print("ERROR: Unknown reconstruction method: " + recon_alg)
        #printTree(tree,0,names=[germanicNameTable[lang] for lang in cognateLangs], field="reconstructed", func="".join)
        cons = "".join(tree.reconstructed)

    if gaps:
        return cons
    else:
        return [c for c in cons if c != '-'] #cons.replace('-','')




#.. File Formats
#.. ------------
#.. 
#.. Pairwise as well as multiple sequence comparison is basically carried out by
#.. reading data from text files and writing the results of the analyses back to
#.. text files. For pairwise and multiple sequence analyses, specific file formats 
#.. are required. See the documentation for the respective classes for details.
#.. 
#.. ``psq``-format
#..     The ``psq``-format is a specific format for text files containing unaligned
#..     sequence pairs. Files in this format should have the extension ``psq``. 
#..     
#..     The first line of a ``psq``-file contains information regarding the dataset.
#..     The sequence pairs are given in triplets, with a sequence identifier in the
#..     first line of a triplet (containing the meaning, or orthographical
#..     information) and the two sequences in the second and third line, whereas
#..     the first column of each sequence line contains the name of the taxon and
#..     the second column the sequence in IPA format. All triplets are divided by
#..     one empty line. As an example, consider the file ``test.psq``::
#.. 
#..         Harry Potter Testset
#..         Woldemort in German and Russian
#..         German  waldemar
#..         Russian vladimir
#.. 
#..         Woldemort in English and Russian
#..         English woldemort
#..         Russian vladimir
#.. 
#..         Woldemort in English and German
#..         English woldemort
#..         German  waldemar
#.. 
#.. ``psa``-format
#..     The ``psa``-format is a specific format for text files containing unaligned
#..     sequence pairs. Files in this format should have the extension ``psq``. 
#..     
#..     The first line of a ``psa``-file contains information regarding the
#..     dataset.  The sequence pairs are given in quadruplets, with a sequence
#..     identifier in the first line of a quadruplet (containing the meaning, or
#..     orthographical information) and the aligned sequences in the second and
#..     third line, whith the name of the taxon in the first column and all aligned
#..     segments in the following columns, separated by tabstops. The fourth line
#..     contains a float indicating the similarity score of the sequences.  All
#..     quadruplets are divided by one empty line. As an example, consider the file
#..     ``test.psa``::
#.. 
#..         Harry Potter Testset
#..         Woldemort in German and Russian
#..         German.    w    a    l    -    d    e    m    a    r
#..         Russian    v    -    l    a    d    i    m    i    r
#..         41.0
#..         
#..         Woldemort in English and Russian
#..         English    w    o    l    -    d    e    m    o    r    t
#..         Russian    v    -    l    a    d    i    m    i    r    -
#..         34.0
#..         
#..         Woldemort in English and German
#..         English    w    o    l    d    e    m    o    r    t
#..         German.    w    a    l    d    e    m    a    r    -
#..         56.0
#.. 
#.. ``msq``-format
#..     The ``msq``-format is a specific format for text files containing unaligned
#..     sequences. Files in this format should have the extension ``msq``. The
#..     first line of an ``msq``-file contains information regarding the dataset.
#..     The second line contains information regarding the sequence (meaning,
#..     identifier), and the following lines contain the name of the taxa in the
#..     first column and the sequences in IPA format in the second column,
#..     separated by a tabstop. As an example, consider the file ``test.msq``::
#.. 
#..         Harry Potter Testset
#..         Woldemort (in different languages)
#..         German  waldemar
#..         English woldemort
#..         Russian vladimir
#.. 
#.. 
#.. ``msa``-format
#..     The ``msa``-format is a specific format for text files containing already
#..     aligned sequence pairs. Files in this format should have the extension
#..     ``msa``. 
#..     
#..     The first line of a ``msa``-file contains information regarding the
#..     dataset. The second line contains information regarding the sequence (its
#..     meaning, the protoform corresponding to the cognate set, etc.). The aligned
#..     sequences are given in the following lines, whereas the taxa are given in
#..     the first column and the aligned segments in the following columns.
#..     Additionally, there may be a specific line indicating the presence of swaps
#..     and a specific line indicating highly consistent sites (local peaks) in the
#..     MSA.  The line for swaps starts with the headword ``SWAPS`` whereas a plus
#..     character (``+``) marks the beginning of a swapped region, the dash
#..     character (``-``) its center and another plus character the end. All sites
#..     which are not affected by swaps contain a dot. The line for local peaks
#..     starts with the headword ``LOCAL``. All sites which are highly consistent
#..     are marked with an asterisk (``*``), all other sites are marked with a dot
#..     (``.``). As an example, consider the file ``test.msa``::
#.. 
#..         Harry Potter Testset
#..         Woldemort (in different languages)
#..         English     w    o    l    -    d    e    m    o    r    t
#..         German.     w    a    l    -    d    e    m    a    r    -
#..         Russian     v    -    l    a    d    i    m    i    r    -
#..         SWAPS..     .    +    -    +    .    .    .    .    .    .
#..         LOCAL..     *    *    *    .    *    *    *    *    *    .


