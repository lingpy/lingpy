# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-13 18:31
# modified : 2013-03-15 10:48
"""
Evaluation methods for automatic cognate detection.
"""

__author__="Johann-Mattis List"
__date__="2013-03-15"

from ..check.messages import FileWriteMessage

def bcubes(
        lex,
        gold='cogid',
        test='lexstatid',
        loans=True,
        pprint=True
        ):
    """
    Compute B-Cubed scores for test and reference datasets.

    Parameters
    ----------
    lex : :py:class:`lingpy.compare.lexstat.LexStat`
        The :py:class:`~lingpy.compare.lexstat.LexStat` class used for the
        computation. It should have two columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    loans : bool (default=True)
        If set to c{False}, loans (indicated by negative IDs in the gold
        standard) will be treated as separate cognates, otherwise, loans will
        be treated as cognates.
    pprint : bool (default=True)
        Print out the results

    Returns
    -------
    t : tuple
        A tuple consisting of the precision, the recall, and the harmonic mean
        (F-scores).

    Notes
    -----
    B-Cubed scores were first described by :evobib:`Bagga1998` as part of an
    algorithm. Later on, :evobib:`Amigo2009` showed that they can also used as
    to compare cluster decisions. :evobib:`Hauer2011` applied the B-Cubed
    scores first to the task of automatic cognate detection.
    
    See also
    --------
    diff
    pairs
    """

    # get the etymdicts
    etdG = []
    for key,line in lex.get_etymdict(ref=gold,loans=loans).items():
        etdG += [[]]
        for value in [x for x in line if x != 0]:
            etdG[-1] += value
    etdT = []
    for key,line in lex.get_etymdict(ref=test,loans=loans).items():
        etdT += [[]]
        for value in [x for x in line if x != 0]:
            etdT[-1] += value
    
    # b-cubed recall
    bcr = []

    # start comparing from gold perspective
    for lineG in etdG:
        
        # get the line of the cluster
        gLen = len(lineG)

        # check for linesize
        if gLen > 1:

            # get cognate-ids in the testset for the line
            lineT = [lex[idx,test] for idx in lineG]
            
            # get the recall
            for idx in lineT:
                bcr += [lineT.count(idx) / gLen]
        
        # otherwise
        else:
            bcr += [1.0]
    
    # b-cubed precision
    bcp = []

    # now compare from the test perspective
    for lineT in etdT:
        
        # get the line of the cluster
        tLen = len(lineT)

        # check for linesize
        if tLen > 1:

            # get cognate-ids in the testset for the line
            lineG = [lex[idx,gold] for idx in lineT]
            
            # get the recall
            for idx in lineG:
                bcp += [lineG.count(idx) / tLen]
        
        # otherwise
        else:
            bcp += [1.0]

    # calculate general scores
    BCP = sum(bcp) / len(bcp)
    BCR = sum(bcr) / len(bcr)
    FSC = 2 * ( ( BCP * BCR ) / ( BCP + BCR ) )

    # print the results if this option is chosen
    if pprint:
        print('*****************************')
        print('* B-Cubed-Scores            *')
        print('* ------------------------- *')
        print('* B-Cubed-Precision: {0:.4f} *'.format(BCP))
        print('* B-Cubed-Recall:    {0:.4f} *'.format(BCR))
        print('* B-Cubed-F-Scores:  {0:.4f} *'.format(FSC))
        print('*****************************')

    # return the stuff
    return BCP,BCR,FSC

def pairs(
        lex,
        gold='cogid',
        test='lexstatid',
        loans=True,
        pprint=True
        ):
    """
    Compute pair scores for the evaluation of cognate detection algorithms.
    
    .. , following Bouchard-Côté et al. (2013).
    lex : :py:class:`lingpy.compare.lexstat.LexStat`
        The :py:class:`~lingpy.compare.lexstat.LexStat` class used for the
        computation. It should have two columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    loans : bool (default=True)
        If set to c{False}, loans (indicated by negative IDs in the gold
        standard) will be treated as separate cognates, otherwise, loans will
        be treated as cognates.
    pprint : bool (default=True)
        Print out the results

    Returns
    -------
    t : tuple
        A tuple consisting of the precision, the recall, and the harmonic mean
        (F-scores).
    
    Notes
    -----
    Pair-scores can be computed in different ways, with often different
    results. This variant follows the description by :evobib:`Bouchard-Cote2013`.
    
    See also
    --------
    diff
    bcubes
    """

    # get the etymdicts
    etdG = []
    for key,line in lex.get_etymdict(ref=gold,loans=loans).items():
        etdG += [[]]
        for value in [x for x in line if x != 0]:
            etdG[-1] += value
    etdT = []
    for key,line in lex.get_etymdict(ref=test,loans=loans).items():
        etdT += [[]]
        for value in [x for x in line if x != 0]:
            etdT[-1] += value
    
    # get the pairs for gold and test
    pairsG = []
    for line in etdG:
        for i,a in enumerate(line):
            for j,b in enumerate(line):
                if i < j:
                    pairsG += [tuple(sorted([a,b]))]
    pairsG = set(pairsG)

    pairsT = []
    for line in etdT:
        for i,a in enumerate(line):
            for j,b in enumerate(line):
                if i < j:
                    pairsT += [tuple(sorted([a,b]))]
    pairsT = set(pairsT)
    
    # calculate precision and recall
    pp = len(pairsG.intersection(pairsT)) / len(pairsT)
    pr = len(pairsG.intersection(pairsT)) / len(pairsG)
    fs = 2 * ( pp * pr ) / ( pp + pr )

    # print the results if this option is chosen
    if pprint:
        print('**************************')
        print('* Pair-Scores            *')
        print('* ---------------------- *')
        print('* Pair-Precision: {0:.4f} *'.format(pp))
        print('* Pair-Recall:    {0:.4f} *'.format(pr))
        print('* Pair-F-Scores:  {0:.4f} *'.format(fs))
        print('**************************')
    
    return pp,pr,fs

def diff(
        lex,
        gold='cogid',
        test='lexstatid',
        loans=True,
        pprint=True,
        filename = '',
        tofile = True
        ):
    r"""
    Write differences in classifications on an item-basis to file.

    lex : :py:class:`lingpy.compare.lexstat.LexStat`
        The :py:class:`~lingpy.compare.lexstat.LexStat` class used for the
        computation. It should have two columns indicating cognate IDs.
    gold : str (default='cogid')
        The name of the column containing the gold standard cognate
        assignments.
    test : str (default='lexstatid')
        The name of the column containing the automatically implemented cognate
        assignments.
    loans : bool (default=True)
        If set to c{False}, loans (indicated by negative IDs in the gold
        standard) will be treated as separate cognates, otherwise, loans will
        be treated as cognates.
    pprint : bool (default=True)
        Print out the results
    filename : str (default='')
        Name of the output file. If not specified, it is identical with the
        name of the :py:class:`~lingpy.compare.lexstat.LexStat`, but with the
        extension ``diff``.
    tofile : bool (default=True)
        If set to c{False}, no data will be written to file, but instead, the
        data will be returned.

    Returns
    -------
    t : tuple
        A nested tuple consisting of two further tuples. The first 
        containing precision, recall, and harmonic mean
        (F-scores), the second containing the same values for the pair-scores.

    Notes
    -----
    If the **tofile** option is chosen, the results are written to a specific
    file with the extension ``diff``. This file contains all cognate sets in
    which there are differences between gold standard and test sets. It also
    gives detailed information regarding false positives, false negatives, and
    the words involved in these wrong decisions.

    .. This function also calculates the "transformation" score. This score is
    .. based on the calculation of steps that are needed to transform one cluster
    .. for one set of meanings into the other. Ideally, if there are *n* different
    .. cognate sets covering one gloss in the gold standard, the minimal length of 
    .. a mapping to convert the *m* cognate sets of the test set into the gold standard
    .. is *n*. In this case, both gold standard and test set are identical.
    .. However, if gold standard and test set differ, the number of mappings
    .. necessarily exceeds *m* and *n*. Based on this, the transformation
    .. precision is defined as :math:`\frac{m}{M}`, where *m* is the number of
    .. distinct clusters in the test set and *M* is the length of the mapping.
    .. Accordingly, the recall is defined as :math:`\frac{n}{M}`, where *n* is the
    .. number of clusters in the gold standard.

    .. Note that if precision is lower than 1.0, this means there are false
    .. positive decisions in the test set. Accordingly, a recall lower than 1.0
    .. indicates that there are false negative decisions in the test set.
    .. The drawback of this score is that it is not sensitive regarding the
    .. distinct number of decisions in which gold standard and test set differ, so
    .. the recall can be very low although most of the words have been grouped
    .. accurately. The advantage is that it can be directly interpreted in terms
    .. of 'false positive/false negative' decisions.

    See also
    --------
    bcubes
    pairs
    """
    # check for filename
    if not filename:
        filename = lex.filename
    else:
        pass

    if loans:
        loan = lambda x: abs(x)
    else:
        loan = lambda x: x

    # open file
    if tofile:
        f = open(filename+'.diff','w')

    # get a formatter for language names
    lform = '{0:'+str(max([len(l) for l in lex.cols]))+'}'
    
    preT,recT = [],[]
    preB,recB = [],[]
    preP,recP = [],[]

    # iterate over all concepts
    for concept in lex.concepts:
        idxs = lex.get_list(row=concept,flat=True)
        cogsG = lex.get_list(row=concept,entry=gold,flat=True)
        cogsT = lex.get_list(row=concept,entry=test,flat=True)

        # compare cogs and test
        # get the basic index for all seqs
        bidx = [i+1 for i in range(len(idxs))]
        
        # translate to cogs
        tmp = {}
        for a,b in zip(cogsG,bidx):
            if loan(a) in tmp:
                pass
            else:
                tmp[loan(a)] = b
        cogsG = [tmp[loan(i)] for i in cogsG]

        # translate to test
        tmp = {}
        for a,b in zip(cogsT,bidx):
            if loan(a) in tmp:
                pass
            else:
                tmp[loan(a)] = b
        cogsT = [tmp[loan(i)] for i in cogsT]
        
        if cogsG != cogsT:

            # calculate the transformation distance of the sets
            tramGT = len(set(zip(cogsG,cogsT)))
            tramG  = len(set(cogsG))
            tramT  = len(set(cogsT))
            preT += [tramT/tramGT]
            recT += [tramG/tramGT]

            # calculate the bcubed precision for the sets
            tmp = {}
            for x,y in zip(cogsT,cogsG):
                if x in tmp:
                    tmp[x] += [y]
                else:
                    tmp[x] = [y]
            bcp = 0.0
            for x in tmp:
                for y in tmp[x]:
                    bcp += tmp[x].count(y) / len(tmp[x])
            preB += [bcp / len(idxs)]
            
            # calculate b-cubed recall
            tmp = {}
            for x,y in zip(cogsG,cogsT):
                if x in tmp:
                    tmp[x] += [y]
                else:
                    tmp[x] = [y]
            bcr = 0.0
            for x in tmp:
                for y in tmp[x]:
                    bcr += tmp[x].count(y) / len(tmp[x])
            recB += [bcr / len(idxs)]

            # calculate pair precision
            tmp = {}
            for x,y in zip(cogsG,idxs):
                if x in tmp:
                    tmp[x] += [y]
                else:
                    tmp[x] = [y]
            pairsG = []
            for x in tmp:
                for i,yA in enumerate(tmp[x]):
                    for j,yB in enumerate(tmp[x]):
                        if i < j:
                            pairsG += [tuple(sorted([yA,yB]))]
            pairsG = set(pairsG)
            
            tmp = {}
            for x,y in zip(cogsT,idxs):
                if x in tmp:
                    tmp[x] += [y]
                else:
                    tmp[x] = [y]
            pairsT = []
            for x in tmp:
                for i,yA in enumerate(tmp[x]):
                    for j,yB in enumerate(tmp[x]):
                        if i < j:
                            pairsT += [tuple(sorted([yA,yB]))]
            pairsT = set(pairsT)
            
            # append stuff
            if pairsT:
                preP += [len(pairsT.intersection(pairsG)) / len(pairsT) ]
            else:
                preP += [1.0]
            if pairsG:
                recP += [len(pairsT.intersection(pairsG)) / len(pairsG) ]
            else:
                recP += [1.0]

            if preP[-1] == 1.0:
                fp = "no"
            else:
                fp = "yes"
            if recP[-1] == 1.0:
                fn = "no"
            else:
                fn = "yes"
            
            if tofile:
                f.write("Concept: {0}, False Positives: {1}, False Negatives: {2}\n".format(
                    concept,
                    fp,
                    fn
                    ))

            # get the words
            words = [lex[i,'ipa'] for i in idxs]
            langs = [lex[i,'taxa'] for i in idxs]

            # get a word-formater
            wform = '{0:'+str(max([len(w) for w in words]))+'}'

            # write differences to file
            if tofile:
                for word,lang,cG,cT in sorted(
                        zip(
                            words,
                            langs,
                            cogsG,
                            cogsT
                            ),
                        key=lambda x:(x[2],x[3])
                        ):
                    f.write('{0}\t{1}\t{2:4}\t{3:4}\n'.format(
                        lform.format(lang),
                        wform.format(word),
                        cG,
                        cT
                        ))
                f.write('#\n')
        else:
            preT += [1.0]
            recT += [1.0]
            preB += [1.0]
            recB += [1.0]
            preP += [1.0]
            recP += [1.0]

    hp = sum(preT) / len(preT)
    hr = sum(recT) / len(recT)
    hf = 2 * ( hp * hr ) / ( hp + hr )
    bp = sum(preB) / len(preB)
    br = sum(recB) / len(recB)
    bf = 2 * ( bp * br ) / ( bp + br )
    pp = sum(preP) / len(preP)
    pr = sum(recP) / len(recP)
    pf = 2 * ( pp * pr ) / ( pp + pr )

    if pprint:
        #print('**************************')
        #print('* Transformation-Scores  *')
        #print('* ---------------------- *')
        #print('* Tram-Precision: {0:.4f} *'.format(hp))
        #print('* Tram-Recall:    {0:.4f} *'.format(hr))
        #print('* Tram-F-Scores:  {0:.4f} *'.format(hf))
        #print('**************************')
        #print('')
        print('**************************')
        print('* B-Cubed-Scores         *')
        print('* ---------------------- *')
        print('* B-C.-Precision: {0:.4f} *'.format(bp))
        print('* B-C.-Recall:    {0:.4f} *'.format(br))
        print('* B-C.-F-Scores:  {0:.4f} *'.format(bf))
        print('**************************')
        print('')
        print('**************************')
        print('* Pair-Scores            *')
        print('* ---------------------- *')
        print('* Pair-Precision: {0:.4f} *'.format(pp))
        print('* Pair-Recall:    {0:.4f} *'.format(pr))
        print('* Pair-F-Scores:  {0:.4f} *'.format(pf))
        print('**************************')
    
    if tofile:
        f.write('B-Cubed Scores:\n')
        f.write('Precision: {0:.4f}\n'.format(bp))
        f.write('Recall:    {0:.4f}\n'.format(br))
        f.write('F-Score:   {0:.4f}\n'.format(bf))
        f.write('#\n')
        f.write('Pair Scores:\n')
        f.write('Precision: {0:.4f}\n'.format(pp))
        f.write('Recall:    {0:.4f}\n'.format(pr))
        f.write('F-Score:   {0:.4f}\n'.format(pf))
        #f.write('#\n')
        #f.write('Hamming Scores:\n')
        #f.write('Precision: {0:.4f}\n'.format(hp))
        #f.write('Recall:    {0:.4f}\n'.format(hr))
        #f.write('F-Score:   {0:.4f}\n'.format(hf))
        #f.write('#\n')
        f.close()
        FileWriteMessage(filename,'diff').message('written')

    else:
        return ((bp,br,bf),(pp,pr,pf))
    
