# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-07 20:07
# modified : 2013-03-07 21:23

"""
Basic module for pairwise and multiple sequence comparison.

The module consists of four classes which deal with pairwise and multiple
sequence comparison from the *sequence* and the *alignment* perspective. The
sequence perspective deals with unaligned sequences. The *alignment*
perspective deals with aligned sequences.

File Formats
------------

Pairwise as well as multiple sequence comparison is basically carried out by
reading data from text files and writing the results of the analyses back to
text files. For pairwise and multiple sequence analyses, specific file formats 
are required. See the documentation for the respective classes for details.

``psq``-format
    The ``psq``-format is a specific format for text files containing unaligned
    sequence pairs. Files in this format should have the extension ``psq``. 
    
    The first line of a ``psq``-file contains information regarding the dataset.
    The sequence pairs are given in triplets, with a sequence identifier in the
    first line of a triplet (containing the meaning, or orthographical
    information) and the two sequences in the second and third line, whereas
    the first column of each sequence line contains the name of the taxon and
    the second column the sequence in IPA format. All triplets are divided by
    one empty line. As an example, consider the file ``test.psq``::

        Harry Potter Testset
        Woldemort in German and Russian
        German  waldemar
        Russian vladimir

        Woldemort in English and Russian
        English woldemort
        Russian vladimir

        Woldemort in English and German
        English woldemort
        German  waldemar

``psa``-format
    The ``psa``-format is a specific format for text files containing unaligned
    sequence pairs. Files in this format should have the extension ``psq``. 
    
    The first line of a ``psa``-file contains information regarding the
    dataset.  The sequence pairs are given in quadruplets, with a sequence
    identifier in the first line of a quadruplet (containing the meaning, or
    orthographical information) and the aligned sequences in the second and
    third line, whith the name of the taxon in the first column and all aligned
    segments in the following columns, separated by tabstops. The fourth line
    contains a float indicating the similarity score of the sequences.  All
    quadruplets are divided by one empty line. As an example, consider the file
    ``test.psa``::

        Harry Potter Testset
        Woldemort in German and Russian
        German.    w    a    l    -    d    e    m    a    r
        Russian    v    -    l    a    d    i    m    i    r
        41.0
        
        Woldemort in English and Russian
        English    w    o    l    -    d    e    m    o    r    t
        Russian    v    -    l    a    d    i    m    i    r    -
        34.0
        
        Woldemort in English and German
        English    w    o    l    d    e    m    o    r    t
        German.    w    a    l    d    e    m    a    r    -
        56.0

``msq``-format
    The ``msq``-format is a specific format for text files containing unaligned
    sequences. Files in this format should have the extension ``msq``. The
    first line of an ``msq``-file contains information regarding the dataset.
    The second line contains information regarding the sequence (meaning,
    identifier), and the following lines contain the name of the taxa in the
    first column and the sequences in IPA format in the second column,
    separated by a tabstop. As an example, consider the file ``test.msq``::

        Harry Potter Testset
        Woldemort (in different languages)
        German  waldemar
        English woldemort
        Russian vladimir


``msa``-format
    The ``msa``-format is a specific format for text files containing already
    aligned sequence pairs. Files in this format should have the extension
    ``msa``. 
    
    The first line of a ``msa``-file contains information regarding the
    dataset. The second line contains information regarding the sequence (its
    meaning, the protoform corresponding to the cognate set, etc.). The aligned
    sequences are given in the following lines, whereas the taxa are given in
    the first column and the aligned segments in the following columns.
    Additionally, there may be a specific line indicating the presence of swaps
    and a specific line indicating highly consistent sites (local peaks) in the
    MSA.  The line for swaps starts with the headword ``SWAPS`` whereas a plus
    character (``+``) marks the beginning of a swapped region, the dash
    character (``-``) its center and another plus character the end. All sites
    which are not affected by swaps contain a dot. The line for local peaks
    starts with the headword ``LOCAL``. All sites which are highly consistent
    are marked with an asterisk (``*``), all other sites are marked with a dot
    (``.``). As an example, consider the file ``test.msa``::

        Harry Potter Testset
        Woldemort (in different languages)
        English     w    o    l    -    d    e    m    o    r    t
        German.     w    a    l    -    d    e    m    a    r    -
        Russian     v    -    l    a    d    i    m    i    r    -
        SWAPS..     .    +    -    +    .    .    .    .    .    .
        LOCAL..     *    *    *    .    *    *    *    *    *    .

"""

__author__="Johann-Mattis List"
__date__="2013-03-07"

import numpy as np
import re
from datetime import datetime

from ..data import *
from ..basic import Wordlist
from ..sequence.sound_classes import *
from .multiple import Multiple
from .pairwise import Pairwise



class _Multiple(Multiple):
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
    In order to read in data from text files, two different file formats can be
    used along with this class:

    *msq-format*
        The ``msq``-format is a specific format for text files containing unaligned
        sequences. Files in this format should have the extension ``msq``. The
        first line of an ``msq``-file contains information regarding the dataset.
        The second line contains information regarding the sequence (meaning,
        identifier), and the following lines contain the name of the taxa in the
        first column and the sequences in IPA format in the second column,
        separated by a tabstop. As an example, consider the file ``test.msq``::
    
            Harry Potter Testset
            Woldemort (in different languages)
            German  waldemar
            English woldemort
            Russian vladimir
    
    *msa-format*
        The ``msa``-format is a specific format for text files containing already
        aligned sequence pairs. Files in this format should have the extension
        ``msa``. 
        
        The first line of a ``msa``-file contains information regarding the
        dataset. The second line contains information regarding the sequence (its
        meaning, the protoform corresponding to the cognate set, etc.). The aligned
        sequences are given in the following lines, whereas the taxa are given in
        the first column and the aligned segments in the following columns.
        Additionally, there may be a specific line indicating the presence of swaps
        and a specific line indicating highly consistent sites (local peaks) in the
        MSA.  The line for swaps starts with the headword ``SWAPS`` whereas a plus
        character (``+``) marks the beginning of a swapped region, the dash
        character (``-``) its center and another plus character the end. All sites
        which are not affected by swaps contain a dot. The line for local peaks
        starts with the headword ``LOCAL``. All sites which are highly consistent
        are marked with an asterisk (``*``), all other sites are marked with a dot
        (``.``). As an example, consider the file ``test.msa``::
    
            Harry Potter Testset
            Woldemort (in different languages)
            English     w    o    l    -    d    e    m    o    r    t
            German.     w    a    l    -    d    e    m    a    r    -
            Russian     v    -    l    a    d    i    m    i    r    -
            SWAPS..     .    +    -    +    .    .    .    .    .    .
            LOCAL..     *    *    *    .    *    *    *    *    *    .

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
            raw_data = open(infile+'.msq','r')
        except IOError:
            raw_data = open(infile,'r')
        
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
            raw_data = open(infile+'.msa','r')
        except:
            raw_data = open(infile,'r')

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
            model = None
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

        self.classes = []
        
        if not model:
            self.model = sca
        else:
            self.model = model

        # redefine the sequences of the Multiple class
        class_strings = [tokens2class(seq.split('.'),self.model)
                for seq in self.seqs]
        
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

        out = open(outfile,'w')

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
            out.write('# Created: {0}\n'.format(datetime.today()))
        except:
            pass
        out.close()

class _Pairwise(Pairwise):
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
    used along with this class:

    *psq-format*
        The ``psq``-format is a specific format for text files containing unaligned
        sequence pairs. Files in this format should have the extension ``psq``. 
        
        The first line of a ``psq``-file contains information regarding the dataset.
        The sequence pairs are given in triplets, with a sequence identifier in the
        first line of a triplet (containing the meaning, or orthographical
        information) and the two sequences in the second and third line, whereas
        the first column of each sequence line contains the name of the taxon and
        the second column the sequence in IPA format. All triplets are divided by
        one empty line. As an example, consider the file ``test.psq``::
    
            Harry Potter Testset
            Woldemort in German and Russian
            German  waldemar
            Russian vladimir
    
            Woldemort in English and Russian
            English woldemort
            Russian vladimir
    
            Woldemort in English and German
            English woldemort
            German  waldemar
    
    *psa-format*
        The ``psa``-format is a specific format for text files containing
        already aligned sequence pairs. Files in this format should have the
        extension ``psq``. 
        
        The first line of a ``psa``-file contains information regarding the
        dataset.  The sequence pairs are given in triplets, with a sequence
        identifier in the first line of a triplet (containing the meaning, or
        orthographical information) and the aligned sequences in the second and
        third line, whith the name of the taxon in the first column and all aligned
        segments in the following columns, separated by tabstops. All
        triplets are divided by one empty line. As an example, consider the file
        ``test.psa``::
    
            Harry Potter Testset
            Woldemort in German and Russian
            German.    w    a    l    -    d    e    m    a    r
            Russian    v    -    l    a    d    i    m    i    r
            
            Woldemort in English and Russian
            English    w    o    l    -    d    e    m    o    r    t
            Russian    v    -    l    a    d    i    m    i    r    -
            
            Woldemort in English and German
            English    w    o    l    d    e    m    o    r    t
            German.    w    a    l    d    e    m    a    r    -

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
            raw_data = open(infile+'.psa','r')
        except:
            raw_data = open(infile,'r')

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
            raw_data = open(infile+'.psq','r')
        except:
            raw_data = open(infile,'r')

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
            tmp = open(outfile)
            tmp.close()
            outfile = filename + '_out.' + fileformat
        except:
            pass

        # open output file
        out = open(outfile,'w')

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

class _SCA(Wordlist):
    """
    Class handles Wordlists for the purpose of alignment analyses.
    """
    def __init__(
            self,
            infile,
            row = 'concept',
            col = 'doculect',
            conf = '',
            **keywords
            ):
        
        # initialize the wordlist
        Wordlist.__init__(self,infile,row,col,conf)
        
        # check for cognate-id or alignment-id in header
        if 'cogid' in self.header or 'almid' in self.header:
            self.etd = self.get_etymdict(loans=True)

        # check for strings
        if 'tokens' in self.header:
            stridx = self.header['tokens']
        elif 'orthoparse' in self.header:
            stridx = self.header['orthoparse']
        elif 'ipa' in self.header:
            stridx = self.header['ipa']
        else:
            print("[i] No valid source for strings could be found")
            return 

        # create the alignments by assembling the ids of all sequences
        if 'msa' not in self._meta:
            self._meta['msa'] = {}
            for key,value in self.etd.items():
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
                    if 'concept' in self.header:
                        concept = self[seqids[0],'concept']
                        d['seq_id'] = 'CogId:{0} ({1})'.format(key,concept)
                    else:
                        d['seq_id'] = 'CogId:{0}'.format(key)
                    
                    # set up the data
                    for seq in seqids:
                        taxon = self[seq,'taxa']
                        string = self[seq][stridx]
                        
                        d['taxa'] += [self[seq,'taxa']]
                        d['seqs'] += [self[seq][stridx]]
                    self._meta['msa'][key] = d

    def align(
            self,
            method = 'progressive',
            iteration = False,
            swap_check = False,
            output = False,
            model = None,
            mode = 'global',
            modes = [
                ('global',-2,0.5),
                ('local',-1,0.5),
                ],
            gop = -3,
            scale = 0.5,
            factor = 0.3,
            tree_calc = 'neighbor',
            gap_weight = 0.5,
            restricted_chars = 'T_',
            classes = True,
            sonar = True,
            scorer = {},
            verbose = True,
            plot = False
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
        for key,value in sorted(self.msa.items(),key=lambda x:x[0]):
            if verbose: print("[i] Analyzing cognate set number {0}.".format(key))

            m = SCA(value)
            if method == 'progressive':
                m.prog_align(
                        model,
                        mode,
                        gop,
                        scale,
                        factor,
                        tree_calc,
                        gap_weight,
                        restricted_chars,
                        classes,
                        sonar,
                        scorer
                        )
            elif method == 'library':
                m.lib_align(
                        model,
                        mode,
                        modes,
                        scale,
                        factor,
                        tree_calc,
                        gap_weight,
                        restricted_chars,
                        classes,
                        sonar,
                        scorer
                        )

            if iteration:
                m.iterate_clusters(0.5)
                m.iterate_orphans()
                m.iterate_similar_gap_sites()

            if swap_check:
                m.swap_check()

            if output:
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
                #self._meta['msa'][key] = {}
            #print(m._alm_matrix)
            self._meta['msa'][key]['alignment'] = m.alm_matrix
                #if plot:
                #    alm2html(
                #            '{0}-msa/{1}-{2}'.format(
                #                self.filename,
                #                self.dataset,
                #                key
                #                )
                #            )

                    
    def __len__(self):
        return len(self.msa)

class SCA(object):
    """
    Basic class for handling various kinds of alignment analyses.

    Parameters
    ----------
    infile : { str }
        The name of the input file which can be provided in PSQ, PSQ, MSQ, or
        MSA format.

    """
    
    def __init__(
            self, 
            infile, 
            **keywords
            ):

        # set the defaults
        defaults = {
                'comment'      : '#',
                "diacritics"   : None,
                "vowels"       : None,
                "tones"        : None,
                "combiners"    : '\u0361\u035c',
                "breaks"       : '.-',
                "stress"       : "ˈˌ'",
                "merge_vowels" : True
                }

        # check for keywords
        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]
        
        # check for datatype
        if type(infile) == dict:
            parent = _Multiple(infile,**keywords)
            case = 1
        else:
            if infile[-4:] in ['.msa','.msq']:
                parent = _Multiple(infile,**keywords)
                case = 1
            elif infile[-4:] in ['.psq','psa']:
                parent = _Pairwise(infile,**keywords)
                case = 2
            elif infile[-4:] in ['.csv','.tsv']:
                parent = _SCA(infile,**keywords)
                case = 3
        
        if case in [1,2]:
            for key in dir(parent):
                if key not in [
                        '__class__',
                        ]:
                    try:
                        setattr(self,key,parent.__getattribute__(key))
                    except:
                        pass
        else:
            for key in dir(parent):
                if not key.startswith('__'):
                    try:
                        setattr(self,key,parent.__getattribute__(key))
                    except:
                        pass

        self._parent = parent
    
    def __len__(self):
        return self._parent.__len__()
    def __str__(self):
        return self._parent.__str__()
    def __repr__(self):
        return self.parent.__repr__()
    def __eq__(self,other):
        return self._parent.__eq__(other)
    def __getitem__(self,x):
        return self._parent.__getitem__(x)
    def __getattr__(self,x):
        return self._parent.__getattr__(x)

        


