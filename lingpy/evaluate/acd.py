# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-03-13 18:31
# modified : 2013-03-13 18:31
"""
Evaluation methods for automatic cognate detection.
"""

__author__="Johann-Mattis List"
__date__="2013-03-13"

def bcubes(
        lex,
        gold='cogid',
        test='lexstatid',
        loans=True,
        pprint=True
        ):
    """
    Compute bcubed-scores for test and reference datasets.
    """

    # get the etymdicts
    etdG = [[x[0] for x in v if x != 0] for k,v in 
            lex.get_etymdict(ref=gold,loans=loans).items()]
    etdT = [[x[0] for x in v if x != 0] for k,v in
            lex.get_etymdict(ref=test,loans=loans).items()]
    
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
        print('Precision: {0:.2f}'.format(BCP))
        print('Recall:    {0:.2f}'.format(BCR))
        print('F-Scores:  {0:.2f}'.format(FSC))

    # return the stuff
    return BCP,BCR,FSC




