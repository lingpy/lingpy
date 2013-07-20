# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-02 06:55
# modified : 2013-06-26 17:40
"""
Miscellaneous routines for data conversion.
"""

__author__="Johann-Mattis List"
__date__="2013-06-26"

from ..settings import rcParams
try:
    from ..algorithm.cython import cluster,misc
except:
    from ..algorithm.cython import _cluster as cluster
    from ..algorithm.cython import _misc as misc

def msa2str(msa,wordlist=False):
    """
    Function converts an MSA object into a string.
    """
    
    out = ''

    # if wordlist ist set to True, don't write the header line and put the
    # after comment
    if wordlist: 
        formatter = "{0}\t{1:"+str(max([len(t) for t in msa['taxa']]))+'}'
        out += '# ' +msa['seq_id']+'\n'
        for a,b,c in zip(msa['ID'],msa['taxa'],msa['alignment']):
            out += formatter.format(a,b)+'\t'
            out += '\t'.join(c)+'\n'
        alm_len = len(c)

    elif type(msa) == dict:
        # get formatter
        formatter = '{0:'+str(max([len(t) for t in msa['taxa']]))+'}'
        out += msa['dataset']+'\n'
        out += msa['seq_id']+'\n'
        for a,b in zip(msa['taxa'],msa['alignment']):
            out += formatter.format(a)+'\t'
            out += '\t'.join(b)+'\n'
        alm_len = len(b)
    else:
        # get formatter
        formatter = '{0:'+str(max([len(t) for t in msa.taxa]))+'}'
        out += msa.dataset+'\n'
        out += msa.seq_id+'\n'
        for a,b in zip(msa.taxa,msa.alm_matrix):
            out += formatter.format(a)+'\t'
            out += '\t'.join(b)+'\n'
        alm_len = len(b)
    
    if 'local' in msa:
        local = msa['local']
    elif hasattr(msa,'local'):
        local = msa.local
    else:
        local = False

    if 'swaps' in msa:
        swaps = msa['swaps']
    elif hasattr(msa,'swaps'):
        swaps = msa.swaps
    else:
        swaps = False
    
    if local:
        if wordlist:
            out += formatter.format(0,"LOCAL")+'\t'
        else:
            out += formatter.format("LOCAL")+'\t'
        tmp = []
        for i in range(alm_len):
            if i in local:
                tmp += ['*']
            else:
                tmp += ['.']
        out += '\t'.join(tmp)+'\n'
    if swaps:
        if wordlist:
            out += formatter.format(0,"SWAPS")+'\t'
        else:
            out += formatter.format("SWAPS")+'\t'
        tmp = alm_len * ['.']
        for swap in swaps:
            a,b,c = swap
            tmp[a] = '+'
            tmp[b] = '-'
            tmp[c] = '+'
        out += '\t'.join(tmp)+'\n'

    return out

def wl2dst(
        wl, # wordlist object
        taxa = "taxa",
        concepts = "concepts",
        cognates = "cogid",
        ):
    """
    Function converts wordlist to distance matrix.
    """
    # check for taxon attribute
    taxa = getattr(wl,taxa)

    # check for concepts
    concepts = getattr(wl,concepts)

    distances = []

    for i,taxA in enumerate(taxa):
        for j,taxB in enumerate(taxa):
            if i < j:
                
                # get the two dictionaries
                dictA = wl.get_dict(col=taxA,entry=cognates)
                dictB = wl.get_dict(col=taxB,entry=cognates)
    
                # count amount of shared concepts
                shared = 0
                missing = 0
                for concept in concepts:
                    try:
                        if [k for k in dictA[concept] if k in dictB[concept]]:
                            shared += 1
                        else:
                            pass
                    except KeyError:
                        missing += 1
    
                # append to distances, catch ZeroDivisionError
                try:
                    distances += [ 1 - shared / (wl.height-missing)]
                except ZeroDivisionError:
                    print(rcParams['E_zero_division'].format(taxA,taxB))
                    distances += [1.0]
    distances = misc.squareform(distances)
    return distances

def matrix2groups(
        threshold,
        distances,
        taxa
        ):
    """
    Calculate flat cluster of distance matrix
    """
    
    flats = cluster.flat_upgma(
            threshold,
            distances,
            revert=True
            )
            
    groups = [flats[i] for i in range(len(taxa))]
    
    # renumber the groups
    groupset = sorted(set(groups))
    renum = dict([(i,j+1) for i,j in zip(
        groupset,
        range(len(groupset))
        )])
    groups = [renum[i] for i in groups]

    return dict(zip(taxa,['G_{0}'.format(g) for g in groups]))

def scorer2str(
        scorer
        ):
    """
    Convert a scoring function to a string.
    """
    
    # get sorted representation of characters
    chars = sorted(
            scorer.chars2int
            )

    # get the matrix
    matrix = scorer.matrix
    
    out = ''

    # write stuff to string
    for i,charA in enumerate(chars):
        out += charA
        for j,charB in enumerate(chars):
            out += '\t{0:.2f}'.format(scorer[charA,charB])
        out += '\n'

    return out

def newick2matrix(newick,labels):
    """
    Function converts a newick-representation of a string into a tree-matrix.
    """

    pass
            
