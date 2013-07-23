# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-03-07 20:07
# modified : 2013-07-18 10:41

"""
Basic module for pairwise and multiple sequence comparison.

The module consists of four classes which deal with pairwise and multiple
sequence comparison from the *sequence* and the *alignment* perspective. The
sequence perspective deals with unaligned sequences. The *alignment*
perspective deals with aligned sequences.

"""

__author__="Johann-Mattis List"
__date__="2013-07-18"

import numpy as np
import re
import codecs
import os

from ..settings import rcParams
from ..basic.wordlist import Wordlist
from ..sequence.sound_classes import *
from .multiple import Multiple
from .pairwise import Pairwise
try:
    from ..algorithm.cython import misc
except:
    from ..algorithm.cython import _misc as misc

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

    Examples
    --------
    Get the path to a file from the testset.

    >>> from lingpy import *
    >>> seq_file = get_file('test.seq')

    Load the file into the Multiple class.

    >>> mult = Multiple(seq_file)

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
    the MSA-format. 

    This class inherits the methods of the
    :py:class:`~lingpy.align.multiple.Multiple` class.

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
            if infile.endswith('.msa'):
                self._init_msa(infile,**keywords)
            else:
                self._init_seq(infile,**keywords)

    def _init_seq(
            self,
            infile,
            **keywords
            ):
        """
        Load an ``msq``-file.
        """

        # import the data from the input file
        self.infile = infile.split('/')[-1].replace('.msq','')
        data = []

        # catch the data from the input file
        try:
            raw_data = codecs.open(infile+'.msq','r','utf-8')
        except IOError:
            raw_data = codecs.open(infile,'r','utf-8')
        
        for line in raw_data:
            if not line.startswith(self.comment):
                data.append(line.strip())
        
        # set the first parameters
        self.dataset = data[0]
        self.seq_id = data[1]
        
        # delete the first lines of the data, since they are no longer needed
        del data[0]
        del data[0]
        
        # extract the rest of the data
        for i in range(len(data)):
            data[i] = data[i].split('\t')

        # strip all whitespace from the data
        for i in range(len(data)):
            for j in range(len(data[i])):
                data[i][j] = data[i][j].strip()

        # delete all lines with empty sequences
        i = 0
        while i < len(data):
            if data[i][1] == '-':
                del data[i]
            else:
                i += 1

        # check the data for inconsistencies
        try:
            data = np.array(data)
        except:
            print("[!] Error in file {0}.".format(infile))
            length = len(data[0])
            print(length)
            for i in range(len(data)):
                if len(data[i]) != length:
                    print("[!] Line {0} in {1} is erroneously coded.".format(
                        i+3,infile))
            return

        self.taxa = data[:,0]

        Multiple.__init__(
                self,
                data[:,1],
                **keywords
                )

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
    
    def _init_msa(
            self,
            infile,
            **keywords
            ):
        """
        Load an ``msa``-file.
        """
        
        # import the data from the input file
        self.infile = infile.split('/')[-1].replace('.msa','')
        data = []
        try:
            raw_data = codecs.open(infile+'.msa','r','utf-8')
        except:
            raw_data = codecs.open(infile,'r','utf-8')

        for line in raw_data:
            if not line.startswith(self.comment):
                data.append(line.strip())
    
        # set the first parameters
        self.dataset = data[0]
        self.seq_id = data[1]
        
        # delete the first lines of the data, since they are no longer needed
        del data[0]
        del data[0]
        
        # split the data
        for i in range(len(data)):
            data[i] = data[i].split('\t')
        
        # check for local alignments and swaps in the data
        tmp_taxa = [line[0].strip('.') for line in data]
        if 'LOCAL' in tmp_taxa:
            local_idx = tmp_taxa.index('LOCAL')
            local = data[local_idx][1:]
            self.peak_idx = [i for i in range(len(local)) if local[i] == '*']
            del data[local_idx]
            del tmp_taxa[local_idx]
        else:
            self.peak_idx = []

        if 'SWAPS' in tmp_taxa:
            swap_idx = tmp_taxa.index('SWAPS')
            swaps = data[swap_idx][1:]
            self.swap_index = []
            i = 0
            while swaps:
                if swaps[0] == '.' or swaps[0] == '': # temporary line
                    swaps.pop(0)
                    i += 1
                elif swaps[0] == '+':
                    swaps.pop(0)
                    swaps.pop(0)
                    swaps.pop(0)
                    self.swap_index.append((i,i+1,i+2))
                    i += 3
            del data[swap_idx]
            del tmp_taxa[swap_idx]
        else:
            self.swap_index = []

        # get the rest of the data
        self.taxa = tmp_taxa[:]
        
        # damn the unicode issues in Python2.6!
        self.alm_matrix = [line[1:] for line in data]
        for i,line in enumerate(self.alm_matrix):
            for j,cell in enumerate(line):
                self.alm_matrix[i][j] = cell

        # join the aligned sequences, keep track of the original tokenization
        # by joining it with dots. That looks ugly and nasty and should
        # probably be checked!
        self.seqs = [
                re.sub(
                    '\.+$',
                    '',
                    re.sub(
                        '^\.+',
                        '',
                        re.sub(
                            '\.\.+',
                            '.',
                            '.'.join(alm).replace('-','')
                        )
                    )
                ) for alm in self.alm_matrix
                ]

        Multiple.__init__(
                self,
                self.seqs,
                **keywords
                )

    def ipa2cls(
            self,
            **keywords
            #model = rcParams['sca']
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
        #if not model:
        #    self.model = sca
        #else:
        #    self.model = model

        # redefine the sequences of the Multiple class
        class_strings = [tokens2class(
            seq.split('.'),
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
            ):
        """
        Write data to file.

        Parameters
        ----------
        fileformat : { 'psa', 'msa', 'msq' }
            Indicate which data should be written to file. Select between:

            * 'psa' -- output of all pairwise alignments in ``psa``-format,
            * 'msa' -- output of the multiple alignment in ``msa``-format, or
            * 'msq' -- output of the multiple sequences in ``msq``-format.

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
        if not filename:
            filename = self.infile

        # define the outfile and check, whether it already exists
        outfile = filename + '.' + fileformat

        # create a specific format string in order to receive taxa of equal
        # length
        mtax = max([len(t) for t in self.taxa])
        txf = '{0:.<'+str(mtax)+'}'

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
            if hasattr(self,'peak_idx'):
                if self.peak_idx:
                    out.write(txf.format("LOCAL")+'\t')
                    tmp = ['.'] * len(self.alm_matrix[0])
                    for i in self.peak_idx:
                        tmp[i] = '*'
                    out.write('\t'.join(tmp)+'\n')

            if hasattr(self,'swap_index'):
                if self.swap_index:
                    out.write(txf.format('SWAPS')+'\t')
                    tmp = ['.'] * len(self.alm_matrix[0])
                    for i in self.swap_index:
                        tmp[i[0]] = '+'
                        tmp[i[1]] = '-'
                        tmp[i[2]] = '+'
                    out.write('\t'.join(tmp)+'\n')
        try:
            out.write('# Created using LingPy-2.0\n')
            out.write('# Parameters: '+self.params+'\n')
            out.write('# Created: {0}\n'.format(rcParams['timestamp']))
        except:
            pass
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
    used along with this class: the PSQ-format, and the PSA-format.

    This class inherits the methods of the
    :py:class:`~lingpy.align.pairwise.Pairwise` class.

    """
    def __init__(
            self,
            infile,
            **keywords
            ):
        
        # set the defaults
        defaults = {
                'comment':'#',
                "diacritics" : None,
                "vowels":None,
                "tones":None,
                "combiners":'\u0361\u035c',
                "breaks":'.-',
                "stress":"ˈˌ'",
                "merge_vowels" : True
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
                            [unicode(a) for a in almA],
                            [unicode(b) for b in almB],
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
            filename =  None
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
    are specified by the "cognates" keyword.

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
                    if 'concept' in self.header:
                        concept = self[seqids[0],'concept']
                        d['seq_id'] = 'Cognate Set: {0} ("{1}")'.format(key,concept)
                    else:
                        d['seq_id'] = 'Cognate Set: {0}'.format(key)
                    
                    # set up the data
                    for seq in seqids:
                        taxon = self[seq,'taxa']
                        string = self[seq][stridx]
                        d['ID'] += [seq]
                        d['taxa'] += [self[seq,'taxa']]
                        d['seqs'] += [self[seq][stridx]]
                    self._meta['msa'][ref][key] = d

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

        plot : bool (default=False)
            Determine whether MSA should be plotted in HTML.
        """
        kw = dict(
                method = 'progressive',
                iteration = False,
                swap_check = False,
                output = False,
                model = rcParams['sca'],
                mode = rcParams['align_mode'],
                modes = rcParams['align_modes'],
                gop = rcParams['align_gop'],
                scale = rcParams['align_scale'],
                factor = rcParams['align_factor'],
                tree_calc = rcParams['align_tree_calc'],
                gap_weight = rcParams['gap_weight'],
                restricted_chars = rcParams['restricted_chars'],
                classes = rcParams['classes'],
                sonar = rcParams['sonar'],
                scorer = rcParams['scorer'],
                plot = False,
                ref = rcParams['ref']
                )
        kw.update(keywords)

        for key,value in sorted(
                self.msa[kw['ref']].items(),
                key=lambda x:x[0]
                ):
            if rcParams['verbose']: print("[i] Analyzing cognate set number {0}.".format(key))

            m = SCA(
                    value,
                    **kw
                    )
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

            if kw['output']:
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
            self._meta['msa'][kw['ref']][key]['alignment'] = m.alm_matrix
            self._meta['msa'][kw['ref']][key]['_sonority_consensus'] = m._sonority_consensus
                    
    def __len__(self):
        return len(self.msa)

    def get_consensus(
            self,
            tree = False,
            gaps = False,
            taxa = False,
            classes = False,
            consensus = 'consensus',
            counterpart = 'ipa',
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
        taxa : {c{list} bool} (default=False)
            If *tree* is chosen as a parameter, specify the taxa in order of the aligned
            strings.
        classes : c{bool} (default=False)
            Specify whether sound classes shall be used to calculate the consensus.
        model : ~lingpy.data.model.Model
            A sound class model according to which the IPA strings shall be
            converted to sound-class strings.

        """
        # determine defaults
        defaults = dict(
                model = rcParams['sca'],
                gap_scale = 1.0,
                ref = rcParams['ref']
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
        
        # go on with the analysis
        cons_dict = {}
        for cog in self.etd[ref]:
            if cog in self.msa[ref]:
                if rcParams['verbose']: print("[i] Analyzing cognate set number '{0}'...".format(cog))
                
                # temporary solution for sound-class integration
                if classes == True:
                    classes = []
                    weights = prosodic_weights(
                            prosodic_string(
                                self.msa[ref][cog]['_sonority_consensus']
                                )
                            )
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
                            taxa = taxa,
                            weights = weights,
                            **keywords
                            )
                    classes = True
                else:
                    cons = get_consensus(
                            self.msa[ref][cog]['alignment'],
                            classes = classes,
                            tree = tree,
                            gaps = gaps,
                            taxa = taxa,
                            **keywords
                            )
            # if there's no msa for a given cognate set, this set is a
            # singleton
            else:
                cons = self[[k[0] for k in self.etd[ref][cog] if k != 0][0],counterpart]
            
            # add consensus to dictionary
            cons_dict[cog] = cons
        
        
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
        fileformat : {'csv', 'tre','nwk','dst', 'taxa', 'starling', 'paps.nex', 'paps.csv'}
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
                filename = rcParams['filename']
                )
        kw.update(keywords)

        # define two vars for convenience
        ref = kw['ref']
        filename = kw['filename']

        if 'cognates' in kw:
            print(rcParams['W_deprecation'].format('cognates','ref'))
            ref = kw['cognates']

        if ref != rcParams['ref']:
            rcParams['ref'] = ref

        if fileformat not in ['alm']:
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
                cogids = [self[i,ref] for i in indices]
                cogids = sorted(set(cogids))
                for cogid in cogids:
                    if cogid in self.msa[ref]:
                        for i,alm in enumerate(self.msa[ref][cogid]['alignment']):
                            taxon = self.msa[ref][cogid]['taxa'][i]
                            seq = self.msa[ref][cogid]['seqs'][i]
                            cid = concept2id[concept]
                            out += '\t'.join(
                                [
                                    str(cogid),
                                    taxon,
                                    concept,
                                    str(cid),
                                    '\t'.join(alm)
                                    ]
                                )+'\n'
                    else:
                        this_idx = [x for x in self.etd[ref][cogid] if x != 0][0][0]
                        taxon = self[this_idx,'taxon']
                        seq = self[this_idx,'tokens']
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

    Returns
    -------
    cons : c{str}
        A consensus string of the given MSA.
    """
    # set defaults
    defaults = dict(
            model = rcParams['sca'],
            gap_scale = 1.0,
            mode = 'majority',
            gap_score = -10,
            weights = [1 for i in range(len(msa[0]))]
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]
    
    # stores the consensus string
    cons = ''

    # transform the matrix
    if hasattr(msa,'alm_matrix'):
        matrix = misc.transpose(msa.alm_matrix)
    else:
        matrix = misc.transpose(msa)

    # check for classes
    if classes:
        # if classes are passed as array, we use this array as is
        if type(classes) == list:
            pass
        # if classes is a Model-object
        elif hasattr(msa,'ipa2cls'):
            msa.ipa2cls(model)
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
                cons += chars[0]
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

                cons += chars[0]

    # otherwise, we use a bottom-up parsimony approach to determine the best
    # match
    elif tree and not taxa:
        raise ValueError(
                "[!] Without a list of taxa, no consensus string can be calculated"
                )
    elif tree and taxa:
        pass # XXX

    if gaps:
        return cons
    else:
        return cons.replace('-','')




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


