# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-02 06:55
# modified : 2013-04-02 06:55
"""
Miscellaneous routines for data conversion.
"""

__author__="Johann-Mattis List"
__date__="2013-04-02"


try:
    from ..algorithm.cython import cluster,misc
except:
    from ..algorithm.cython import _cluster as cluster
    from ..algorithm.cython import _misc as misc

def msa2str(msa):
    
    out = ''
    if type(msa) == dict:
        # get formatter
        formatter = '{0:'+str(max([len(t) for t in msa['taxa']]))+'}'
        out += msa['dataset']+'\n'
        out += msa['seq_id']+'\n'
        for a,b in zip(msa['taxa'],msa['alignment']):
            out += formatter.format(a)+'\t'
            out += '\t'.join(b)+'\n'
    else:
        # get formatter
        formatter = '{0:'+str(max([len(t) for t in msa['taxa']]))+'}'
        out += msa.dataset+'\n'
        out += msa.seq_id+'\n'
        for a,b in zip(msa.taxa,msa.alm_matrix):
            out += formatter.format(a)+'\t'
            out += '\t'.join(b)+'\n'
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
    
                # append to distances
                distances += [ 1 - shared / (wl.height-missing)]
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


