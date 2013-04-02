# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-04-02 07:01
# modified : 2013-04-02 07:01
"""
Functions for tree calculations and working with Newick files.
"""

__author__="Johann-Mattis List"
__date__="2013-04-02"

import xml.dom.minidom as minidom

#import networkx as nx
from ..thirdparty import cogent as cg
try:
    from ..algorithm.cython import cluster
except:
    from ..algorithm.cython import _cluster as cluster

def xml2nwk(
        infile,
        filename = None
        ):
    """
    Convert xml-based MultiTree format to Newick-format.
    """
    
    # parse the xml-file
    document = {}

    # get the document
    document['document'] = minidom.parse(infile)
    
    # get the hash
    document['hash'] = document['document'].getElementsByTagName('hash')[0]
    
    # get the tree
    document['tree'] = document['hash'].getElementsByTagName('tree')[0]
    
    # get the root
    document['root'] = document['tree'].getElementsByTagName('root')[0]
    

    # now start iteration
    nwk = {0:[]}
    queue = [(document['root'],0)]
    taxa = []
    while queue:
        
        root,idx = queue.pop()
        
        max_idx = max([k for k in nwk if type(k) == int])
    
        try:
            nwk[idx]
        except:
            nwk[idx] = []
    
        # get the children
        children = [c for c in root.childNodes if c.nodeName == 'children']
    
        # get the childs
        childs = [c for c in children[0].childNodes if c.nodeName == 'child'] 
    
        #print("Idx {0} has {1} childs".format(idx,len(childs)))
        
        if childs:
            # iterate over childs
            for i,child in enumerate(childs):
                
                queue += [(child,max_idx+i+1)]
                nwk[idx] += [max_idx+i+1]
                #print("\tAdded {1} to {0}.".format(idx,max_idx+i+1))
    
        else:
            name = [c for c in root.childNodes if c.nodeName == 'pri-name'][0]
            name = name.childNodes[0].data
            
            nwk[idx] = [name]
            taxa.append(name)
    
    # now that a dictionary representation of the tree has been created,
    # convert everything to newick

    # first, create a specific newick-dictionary
    newick = {}

    for i in range(len(nwk)):
        
        #create format-string for children
        children = ['{{x_{0}}}'.format(c) for c in nwk[i]]
    
        # create dictionary to replace previous format string
        if len(children) > 1:
            newick['x_'+str(i)] = '('+','.join(children)+')'
        else:
            newick['x_'+str(i)] = children[0]
    
    # add the taxa
    for taxon in taxa:
        newick['x_'+str(taxon)] = taxon
    
    # create the newick-string
    newick_string = "{x_0};"
    newick_check = newick_string
    
    # start conversion
    i = 0
    while True:
        
        newick_string = newick_string.format(**newick)
        if newick_check == newick_string:
            break
        else:
            newick_check = newick_string
    
    if not filename:
        return newick_string
    else:
        f = open(filename+'.nwk','w')
        f.write(newick_string)
        f.close()
        print("[i] Data has been written to file <{0}.nwk>".format(filename))
        return

def matrix2tree(
        matrix,
        taxa,
        tree_calc = "neighbor",
        distances = True
        ):
    """
    Calculate a tree of a given distance matrix.
    """
    
    if tree_calc == 'upgma':
        algorithm = cluster.upgma
    elif tree_calc == 'neighbor':
        algorithm = cluster.neighbor

    newick = algorithm(matrix,taxa,distances)

    tree = cg.LoadTree(treestring=newick)

    return tree


