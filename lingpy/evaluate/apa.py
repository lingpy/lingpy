# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-02 12:45
# modified : 2013-09-02 13:32
"""
Basic module for the comparison of automatic phonetic alignments.
"""

__author__="Johann-Mattis List"
__date__="2013-09-02"

# external modules
import numpy as np
import codecs

# internal modules
from ..settings import rcParams
from ..align.sca import PSA,MSA
from ..algorithm import misc
from .. import log


class EvalMSA(object):
    """
    Base class for the evaluation of automatic multiple sequence analyses.

    Parameters
    ----------

    gold, test : :py:class:`~lingpy.align.sca.MSA`
        The :py:class:`~lingpy.compare.Multiple` objects which shall be
        compared. The first object should be the gold standard and the second
        object should be the test set. 

    Notes
    -----

    Most of the scores which can be calculated with help of this class are standard
    evaluation scores in evolutionary biology. For a close description on how
    these scores are calculated, see, for example, :evobib:`Thompson1999`,
    :evobib:`List2012`, and :evobib:`Rosenberg2009b`.

    See also
    --------
    lingpy.evaluate.apa.EvalPSA
    """
    def __init__(self, 
            gold, 
            test
            ):
        self.gold = gold
        self.test = test
    
    def _get_c_scores(self):
        """
        Calculate the c-scores.
        """
        
        almsGold = misc.transpose(self.gold.alm_matrix)
        almsTest = misc.transpose(self.test.alm_matrix)

        commons = len([i for i in almsGold if i in almsTest])
        
        self.cp = commons / len(almsTest)
        self.cr = commons / len(almsGold)
        self.c_ = 2 * commons / (len(almsTest) + len(almsGold))
        try:
            self.cf = 2 * self.cp * self.cr / (self.cp + self.cr)
        except ZeroDivisionError:
            self.cf = 0.0

    def c_score(self,mode=1):
        r"""
        Calculate the column (C) score. 

        Parameters
        ----------

        mode : { 1, 2, 3, 4 }
            Indicate, which mode to compute. Select between:

            1. divide the number of common columns in reference and test
               alignment by the total number of columns in the test alignment
               (the traditional C score described in :evobib:`Thompson1999`,
               also known as "precision" score in applications of information
               retrieval),
            
            2. divide the number of common columns in reference and test
               alignment by the total number of columns in the reference
               alignment (also known as "recall" score in applications of
               information retrieval),
            
            3. divide the number of common columns in reference and test
               alignment by the average number of columns in reference and test
               alignment, or

            4. combine the scores of mode ``1`` and mode ``2`` by computing
               their F-score, using the formula :math:`2 * \frac{pr}{p+r}`,
               where *p* is the precision (mode ``1``) and *r* is the recall
               (mode ``2``).

        Returns
        -------
        score : float
            The C score for reference and test alignments.
        
        Notes
        -----
        The different c-

        See also
        --------
        lingpy.evaluate.apa.EvalPSA.c_score

        """
        try:
            self.cp
        except:
            self._get_c_scores()

        if mode == 1:
            return self.cp
        elif mode == 2:
            return self.cr
        elif mode == 3:
            return self.c_
        elif mode == 4:
            return self.cf
        else:
            raise ValueError('The mode you chose is not available.')

    def r_score(
            self
            ):
        """
        Compute the rows (R) score.

        Returns
        -------

        score : float
            The PIR score.

        Notes
        -----
        The R score is the number of identical rows (sequences) in reference and test
        alignment divided by the total number of rows.

        See also
        --------
        lingpy.evaluate.apa.EvalPSA.r_score
        """
        almsGold = self.gold.alm_matrix
        almsTest = self.test.alm_matrix
        
        goods = 0.0
        count = 0.0
        for i in range(len(almsGold)):
            if ''.join(almsGold[i]) == ''.join(almsTest[i]):
                goods += 1.0
            
            count += 1.0

        return goods / count

    def sp_score(self,mode=1):
        """
        Calculate the sum-of-pairs (SP) score.

        Parameters
        ----------

        mode : { 1, 2, 3 }
            Indicate, which mode to compute. Select between:

            1. divide the number of common residue pairs in reference and test
               alignment by the total number of residue pairs in the test
               alignment (the traditional SP score described in
               :evobib:`Thompson1999`, also known as "precision" score in
               applications of information retrieval),
            
            2. divide the number of common residue pairs in reference and test
               alignment by the total number of residue pairs in the reference
               alignment (also known as "recall" score in applications of
               information retrieval),
            
            3. divide the number of common residue pairs in reference and test
               alignment by the average number of residue pairs in reference
               and test alignment.

        Returns
        -------

        score : float
            The SP score for gold standard and test alignments.

        Notes
        -----

        The SP score (see :evobib:`Thompson1999`) is calculated by dividing the number of
        identical residue pairs in reference and test alignment by the total
        number of residue pairs in the reference alignment. 

        See also
        --------
        lingpy.evaluate.apa.EvalPSA.sp_score
        """
        try:
            return self.sp
        except:
            self._pair_scores()
            return self.sp

    def jc_score(self):
        """
        Calculate the Jaccard (JC) score.

        Returns
        -------
        score : float
            The JC score.

        Notes
        -----
        The Jaccard score (see :evobib:`List2012`) is calculated by dividing the size of
        the intersection of residue pairs in reference and test alignment by
        the size of the union of residue pairs in reference and test alignment.

        See also
        --------
        lingpy.test.evaluate.EvalPSA.jc_score

        """
        
        try:
            return self.jc
        except:
            self._pair_scores()
            return self.jc

    def _pair_scores(
            self,
            weights = False
            ):
        """
        Calculate msa alignment scores by calculating the pairwise scores.
        """
        
        if self.gold == self.test:
            self.sp = 1.0 
            self.o1 = 1.0 
            self.o2 = 1.0 
            self.o_ = 1.0 
            self.jc = 1.0 
            self.cg1 = 1.0 
            self.cg2 = 1.0 
            self.cg_ = 1.0 
            self.cgf = 1.0
            self.pip = 1.0
            return 

        # replace all characters by numbers
        almsGold = np.zeros(
                (
                    len(
                        self.gold.alm_matrix
                        ),
                    len(
                        self.gold.alm_matrix[0]
                        )
                    )
                )
        almsTest = np.zeros(
                (
                    len(
                        self.test.alm_matrix
                        ),
                    len(
                        self.test.alm_matrix[0]
                        )
                    )
                )
        
        # select between calculation which is based on an explicit weighting or
        # a calculation which is based on implicit weighting, explicit
        # weighting is done by choosing a specific sound class model and
        # cluster all sequences which are identical, implicit weighting is
        # done otherwise, i.e. identical (pid = 100) sequences are clustered
        # into one sequence in order to avoid getting good scores when there
        # are too many highly identical sequences. 
        # XXX this part of the calculation has never really been testend. I
        # leave it untouched for the moment, since it won't be activated,
        # anyway, but we should come back to this and either follow up the idea
        # or discard the application of weights XXX
        if weights:
            self.gold._set_model(weights)
            self._uniseqs = self.gold.int2ext
        
        else:
            self._uniseqs = {}
            for i,seq in enumerate(self.gold.alm_matrix):
                seq = ''.join(seq).replace('-','')
                try:
                    self._uniseqs[seq] += [i]
                except:
                    self._uniseqs[seq] = [i]

        self.weights = {}
        for key in self._uniseqs.keys():
            vals = self._uniseqs[key]
            l = len(vals)
            for val in vals:
                self.weights[val] = (key,1.0 / l)

        # change residues by assining each residue a unique status in both MSAs
        for key in self._uniseqs:
            k = 1
            vals = self._uniseqs[key]
            tmp = []
            for res in self.gold.alm_matrix[vals[0]]:
                if res == '-':
                    tmp.append(0)
                else:
                    tmp.append(k)
                    k += 1
            almsGold[vals[0]] += np.array(tmp)

        for key in self._uniseqs:
            k = 1
            vals = self._uniseqs[key]
            tmp = []
            for res in self.test.alm_matrix[vals[0]]:
                if res == '-':
                    tmp.append(0)
                else:
                    tmp.append(k)
                    k += 1
            almsTest[vals[0]] += np.array(tmp)

        # start computation by assigning the variables
        crp = 0.0 # common residue pairs
        trp = 0.0 # residue pairs in test alignment
        rrp = 0.0 # residue pairs in reference alignment
        urp = 0.0 # unique residue pairs in test and reference
        gcrp = 0.0 # common residue pairs including gaps
        gtrp = 0.0 # length of test alignment
        grrp = 0.0 # length of reference alignment
        pip = 0.0 # percentage of identical pairs score

        testL = len(almsTest[0])
        goldL = len(almsGold[0])
        
        # start iteration
        for i,almA in enumerate(almsGold):
            for j,almB in enumerate(almsGold):
                if i < j:
                    gold = list(zip(almA,almB))
                    test = list(zip(almsTest[i],almsTest[j]))
                    
                    if self.weights[i][0] != self.weights[j][0]:
                        w = self.weights[i][1] * self.weights[j][1]
                    else:
                        w = 0.0

                    # speed up the stuff when sequences are identical
                    if gold == test:
                        tmp = len([x for x in gold if 0 not in x]) * w
                        crp += tmp 
                        trp += tmp 
                        rrp += tmp 
                        urp += tmp 
                        gcrp += testL * w
                        gtrp += testL * w
                        grrp += goldL * w
                        pip += 1 * w
                    
                    else:
                        if [x for x in gold if x != (0,0)] == [y for y in test \
                                if y != (0,0)]:
                            pip += 1 * w

                        crp += len([x for x in gold if x in test and 0 not in
                            x]) * w
                        trp += len([x for x in test if 0 not in x]) * w
                        rrp += len([x for x in gold if 0 not in x]) * w
                        urp += len(set([x for x in gold+test if 0 not in x])) * w
                        gcrp += len([x for x in gold if x in test]) * w
                        gtrp += testL * w
                        grrp += goldL * w
        
        # calculate the scores
        self.sp = crp / rrp
        self.o1 = self.sp
        self.o2 = crp / trp
        self.o_ = 2 * crp / (rrp + trp)
        self.jc = crp / urp
        self.cg1 = gcrp / grrp # recall
        self.cg2 = gcrp / gtrp # precision
        self.cg_ = 2 * gcrp / (grrp + gtrp)
        self.cgf = 2 * (self.cg1 * self.cg2) / (self.cg1 + self.cg2)
        
        l = len(self._uniseqs)
        self.pip = pip / ((l ** 2 - l) / 2)

    def check_swaps(
            self
            ):
        """
        Check for possibly identical swapped sites.

        Returns
        -------

        swap : { -2, -1, 0, 1, 2 }
            Information regarding the identity of swap decisions is coded by
            integers, whereas

            1 -- indicates that swaps are detected in both gold standard and
              testset, whereas a negative value indicates that the positions
              are not identical,

            2 -- indicates that swap decisions are not identical in gold
              standard and testset, whereas a negative value indicates that
              there is a false positive in the testset, and

            0 -- indicates that there are no swaps in the gold standard and the
              testset.
        """
        try:
            swA = self.gold.swap_index
        except:
            swA = False
        try:
            swB = self.test.swap_index
        except:
            swB = False

        if swA and not swB:
            return 2
        elif not swA and swB:
            return -2
        elif swA and swB:
            if swA == swB:
                return 1
            else:# swA != swB:
                return -1
        elif not swA and not swB:
            return 0

class EvalPSA(object):
    """
    Base class for the evaluation of automatic pairwise sequence analyses.

    Parameters
    ----------
    
    gold, test : :py:class:`lingpy.align.sca.PSA`
        The :py:class:`Pairwise <lingpy.compare.Pairwise>` objects which shall be
        compared. The first object should be the gold standard and the second
        object should be the test set.    
    
    Notes
    -----

    Moste of the scores which can be calculated with help of this class are standard
    evaluation scores in evolutionary biology. For a close description on how
    these scores are calculated, see, for example, :evobib:`Thompson1999`,
    :evobib:`List2012`, and :evobib:`Rosenberg2009b`.

    See also
    --------
    lingpy.evaluate.apa.EvalMSA
    """
    def __init__(self,gold,test):

        self.gold = gold
        self.test = test

    def r_score(
            self,
            mode = 1
            ):
        """
        Compute the percentage of identical rows (PIR) score.

        Parameters
        ----------

        mode : { 1, 2 }
            Select between mode ``1``, where all sequences are compared with
            each other, and mode ``2``, where only whole alignments are
            compared.  

        Returns
        -------

        score : float
            The PIR score.

        Notes
        -----
        The PIR score is the number of identical rows (sequences) in reference and test
        alignment divided by the total number of rows.

        See also
        --------
        lingpy.evaluate.apa.EvalMSA.r_score
        """
        
        score = 0.0
        count = 0.0

        if mode == 1:
            for i,alms in enumerate(self.gold.alignments):
                if self.test.alignments[i][0] == alms[0]:
                    score += 0.5
                if self.test.alignments[i][1] == alms[1]:
                    score += 0.5
                count += 1.0
        elif mode == 2:
            for i,alms in enumerate(self.gold.alignments):
                tmp = 0
                if self.test.alignments[i][0] == alms[0]:
                    tmp = 1
                if self.test.alignments[i][1] == alms[1]:
                    tmp += 1
                if tmp == 2:
                    score += 1.0
                count += 1.0

        return score / count

    def _pairwise_column_scores(
            self
            ):
        """
        Compute the different column scores for pairwise alignments. The method
        returns the precision, the recall score, and the f-score, following the
        proposal of Bergsma and Kondrak (2007), and the column score proposed
        by Thompson et al. (1999).
        """

        # the variables which store the different counts
        
        crp = 0.0 # number of common residue pairs in reference and test alm.
        rrp = 0.0 # number of residue pairs in reference alignment 
        trp = 0.0 # number of residue pairs in test alignment
        urp = 0.0 # number of unique residue pairs in reference and test alm.
        
        gtrp = 0.0 # number of residue pairs (including gaps) in test alm.
        grrp = 0.0 # number of residue pairs (including gaps) in reference alm.
        gcrp = 0.0 # number of common residue pairs (including gaps) in r and t

        self.sps_list = []
        self.cs_list = []

        for i,alms in enumerate(self.gold.alignments):
            zipsA = zip(
                    self.gold.alignments[i][0],
                    self.gold.alignments[i][1]
                    )
            zipsB = zip(
                    self.test.alignments[i][0],
                    self.test.alignments[i][1]
                    )
            
            # replace all residues in reference and test alignment with ids
            pairsGold = []
            j,k = 1,1
            for a,b in zipsA:
                x,y = 0,0
                if a != '-':
                    x = j
                    j += 1
                if b != '-':
                    y = k
                    k += 1
                pairsGold.append((x,y))

            pairsTest = []
            j,k = 1,1
            for a,b in zipsB:
                x,y = 0,0
                if a != '-':
                    x = j
                    j += 1
                if b != '-':
                    y = k
                    k += 1
                pairsTest.append((x,y))

            # calculate the number of residues in crp, rrp, and trp
            commons = len([x for x in pairsTest if x in pairsGold and 0 not in x])
            nogaps = len([x for x in pairsGold if 0 not in x])
            crp += len([x for x in pairsTest if x in pairsGold and 0 not in x])
            rrp += len([x for x in pairsGold if 0 not in x])
            trp += len([x for x in pairsTest if 0 not in x])
            urp += len(set([x for x in pairsGold+pairsTest if 0 not in x]))

            grrp += len(pairsGold)
            gtrp += len(pairsTest)
            gcrp += len([x for x in pairsTest if x in pairsGold])
            
            # fill in list with exact scores
            commons = len([x for x in pairsTest if x in pairsGold and 0 not in x])
            nogaps = len([x for x in pairsGold if 0 not in x])
            nogaps2 = nogaps + len([x for x in pairsTest if 0 not in x])
            columns = len([x for x in pairsTest if x in pairsGold])

            if nogaps != 0:
                self.sps_list.append((2 * commons) / (nogaps + nogaps2)) #new
            elif nogaps == 0 and commons == 0:
                self.sps_list.append(1)
            else:
                self.sps_list.append(0)
            self.cs_list.append((2 * columns)/(len(pairsGold)+len(pairsTest)))#new
        
        # calculate the scores        
        self.sop = crp / rrp
        self.jac = crp / urp
        self.o1 = self.sop
        self.o2 = crp / trp
        self.o_ = 2 * crp / (rrp + trp)
        self.precision = gcrp / gtrp
        self.recall = gcrp / grrp
        self.pic = self.precision
        self.fscore = 2 * (self.precision * self.recall) / (self.precision + \
                self.recall)

    def c_score(
            self
            ):
        """
        Calculate column (C) score. 

        Returns
        -------
        score : float
            The C score for reference and test alignments.
        
        Notes
        -----
        The C score, as it is described in :evobib:`Thompson1999`, is calculated by
        dividing the number of columns which are identical in the gold
        standarad and the test alignment by the total number of columns in the
        test alignment.

        See also
        --------
        lingpy.test.evaluate.EvalMSA.c_score

        """
        try:
            return self.pic
        except:
            self._pairwise_column_scores()
            return self.pic

    def sp_score(
            self
            ):
        """
        Calculate the sum-of-pairs (SP) score.

        Returns
        -------

        score : float
            The SP score for reference and test alignments.

        Note
        ----
        The SP score (see :evobib:`Thompson1999`) is calculated by dividing the number of
        identical residue pairs in reference and test alignment by the total
        number of residue pairs in the reference alignment. 
        
        See also
        --------
        lingpy.test.evaluate.EvalMSA.sp_score

        """
        try:
            return self.sop
        except:
            self._pairwise_column_scores()
            return self.sop

    def jc_score(
            self
            ):
        """
        Calculate the Jaccard (JC) score.

        Returns
        -------
        score : float
            The JC score.

        Notes
        -----

        The Jaccard score (see :evobib:`List2012`) is calculated by dividing the size of
        the intersection of residue pairs in reference and test alignment by
        the size of the union of residue pairs in reference and test alignment.
        
        See also
        --------
        lingpy.test.evaluate.EvalMSA.jc_score

        """
        try:
            return self.jac
        except:
            self._pairwise_column_scores()
            return self.jac

    def diff(
            self,
            **keywords
            ):
        """
        Write all differences between two sets to a file.

        Parameters
        ----------

        filename : str (default='eval_psa_diff')
            Default

        """
        # set up the defaults
        defaults = dict(
                filename = self.gold.infile
                )

        for k in defaults:
            if k not in keywords:
                keywords[k] = defaults[k]

        if not keywords['filename'].endswith('.diff'):
            keywords['filename'] = keywords['filename']+'.diff'

        out = codecs.open(keywords['filename'],'w')

        for i,(a,b) in enumerate(zip(self.gold.alignments,self.test.alignments)):
            
            g1,g2,g3 = a
            t1,t2,t3 = b
            maxL = max([len(g1),len(t1)])
            if g1 != t1 or g2 != t2:
                taxA,taxB = self.gold.taxa[i]
                taxlen = max(len(taxA),len(taxB))
                seq_id = self.gold.seq_ids[i]
                out.write('{0}\n{1}\t{2}\n{3}\t{4}\n{5}\n{1}\t{6}\n{3}\t{7}\n\n'.format(
                    seq_id,
                    taxA,
                    '\t'.join(g1),
                    taxB,
                    '\t'.join(g2),
                    '{0}\t{1}'.format(taxlen*' ','\t'.join(['==' for x in range(maxL)])),
                    '\t'.join(t1),
                    '\t'.join(t2),
                    ))
        out.close()
        log.file_written(keywords['filename'])
