#!/usr/bin/env python
# encoding: utf-8

# link_clustering.py
# Jim Bagrow, Yong-Yeol Ahn
# Last Modified: 2010-08-27

# Copyright 2008,2009,2010 James Bagrow, Yong-Yeol Ahn
# 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
changes 2010-08-27:
    * all three output files now contain the same community id numbers
    * comm2nodes and comm2edges both present the cid as the first
    entry of each line.  Previously only comm2nodes did this.
    * implemented weighted version, added '-w' switch
    * expanded help string to explain input and outputs
"""

# this code was taken and modified from
# https://github.com/bagrow/linkcomm/blob/master/python/link_clustering.py
# accessed on 07.10.2013 (JML)

import sys, os
from copy import copy
from operator import itemgetter
from heapq import heappush, heappop
from collections import defaultdict
from itertools import combinations, chain # requires python 2.6+
from optparse import OptionParser

def swap(a,b):
    if a > b:
        return b,a
    return a,b


def Dc(m,n):
    """partition density"""
    try:
        return m*(m-n+1.0)/(n-2.0)/(n-1.0)
    except ZeroDivisionError: # numerator is "strongly zero"
        return 0.0


class HLC:
    def __init__(self,adj,edges):
        self.adj   = adj # node -> set of neighbors
        self.edges = edges # list of edges
        self.Mfactor  = 2.0 / len(edges)
        self.edge2cid = {}
        self.cid2nodes,self.cid2edges = {},{}
        self.initialize_edges() # every edge in its own comm
        self.D = 0.0 # partition density
    
    def initialize_edges(self):
        for cid,edge in enumerate(self.edges):
            edge = swap(*edge) # just in case
            self.edge2cid[edge] = cid
            self.cid2edges[cid] = set([edge])
            self.cid2nodes[cid] = set( edge )
    
    def merge_comms(self,edge1,edge2):
        if not edge1 or not edge2: # We'll get (None, None) at the end of clustering
            return
        cid1,cid2 = self.edge2cid[edge1],self.edge2cid[edge2]
        if cid1 == cid2: # already merged!
            return
        m1,m2 = len(self.cid2edges[cid1]),len(self.cid2edges[cid2])
        n1,n2 = len(self.cid2nodes[cid1]),len(self.cid2nodes[cid2])
        Dc1, Dc2 = Dc(m1,n1), Dc(m2,n2)
        if m2 > m1: # merge smaller into larger
            cid1,cid2 = cid2,cid1
        
        self.cid2edges[cid1] |= self.cid2edges[cid2]
        for e in self.cid2edges[cid2]: # move edges,nodes from cid2 to cid1
            self.cid2nodes[cid1] |= set( e )
            self.edge2cid[e] = cid1
        del self.cid2edges[cid2], self.cid2nodes[cid2]
        
        m,n = len(self.cid2edges[cid1]),len(self.cid2nodes[cid1]) 
        Dc12 = Dc(m,n)
        self.D = self.D + ( Dc12 -Dc1 - Dc2) * self.Mfactor # update partition density
    
    def single_linkage(self, threshold=None, w=None):
        # print "clustering..." excluded by JML
        self.list_D = [(1.0,0.0)] # list of (S_i,D_i) tuples...
        self.best_D = 0.0
        self.best_S = 1.0 # similarity threshold at best_D
        self.best_P = None # best partition, dict: edge -> cid
        
        if w == None: # unweighted
            H = similarities_unweighted( self.adj ) # min-heap ordered by 1-s
        else: 
            H = similarities_weighted( self.adj, w )
        S_prev = -1
        
        # (1.0, (None, None)) takes care of the special case where the last
        # merging gives the maximum partition density (e.g. a single clique). 
        for oms,eij_eik in chain(H, [(1.0, (None, None))] ):
            S = 1-oms # remember, H is a min-heap
            if threshold:
                if S < threshold:
                    break
                
            if S != S_prev: # update list
                if self.D >= self.best_D: # check PREVIOUS merger, because that's
                    self.best_D = self.D  # the end of the tie
                    self.best_S = S
                    self.best_P = copy(self.edge2cid) # slow...
                self.list_D.append( (S,self.D) )
                S_prev = S
            self.merge_comms( *eij_eik )
        
        #self.list_D.append( (0.0,self.list_D[-1][1]) ) # add final val
        if threshold != None:
            return self.edge2cid, self.D
        return self.best_P, self.best_S, self.best_D, self.list_D
    


def similarities_unweighted(adj):
    """Get all the edge similarities. Input dict maps nodes to sets of neighbors.
    Output is a list of decorated edge-pairs, (1-sim,eij,eik), ordered by similarity.
    """
    # print "computing similarities..." excluded by JML
    i_adj = dict( (n,adj[n] | set([n])) for n in adj)  # node -> inclusive neighbors
    min_heap = [] # elements are (1-sim,eij,eik)
    for n in adj: # n is the shared node
        if len(adj[n]) > 1:
            for i,j in combinations(adj[n],2): # all unordered pairs of neighbors
                edge_pair = swap( swap(i,n),swap(j,n) )
                inc_ns_i,inc_ns_j = i_adj[i],i_adj[j] # inclusive neighbors
                S = 1.0 * len(inc_ns_i&inc_ns_j) / len(inc_ns_i|inc_ns_j) # Jacc similarity...
                heappush( min_heap, (1-S,edge_pair) )
    return [ heappop(min_heap) for i in range(len(min_heap)) ] # return ordered edge pairs


def similarities_weighted(adj, ij2wij):
    """Same as similarities_unweighted but using tanimoto coefficient. `adj' is a dict
    mapping nodes to sets of neighbors, ij2wij is a dict mapping an edge (ni,nj) tuple
    to the weight wij of that edge.  
    """
    # print "computing similarities..." excluded by JML
    i_adj = dict( ( n, adj[n]|set([n]) ) for n in adj ) # node -> inclusive neighbors
    
    Aij = copy(ij2wij)
    n2a_sqrd = {}
    for n in adj:
        if adj[n]: # added check whether node has neighbors or not JML
            Aij[n,n] = 1.0*sum( ij2wij[swap(n,i)] for i in adj[n] )/len(adj[n])
            n2a_sqrd[n] = sum( Aij[swap(n,i)]**2 for i in i_adj[n] ) # includes (n,n)!
    
    min_heap = [] # elements are (1-sim,eij,eik)
    for ind,n in enumerate(adj): # n is the shared node
        #print ind, 100.0*ind/len(adj)
        if len(adj[n]) > 1:
            for i,j in combinations(adj[n],2): # all unordered pairs of neighbors
                edge_pair = swap( swap(i,n),swap(j,n) )
                inc_ns_i,inc_ns_j = i_adj[i],i_adj[j] # inclusive neighbors
                
                ai_dot_aj = 1.0*sum( Aij[swap(i,x)]*Aij[swap(j,x)] for x in inc_ns_i&inc_ns_j )
                
                S = ai_dot_aj / (n2a_sqrd[i]+n2a_sqrd[j]-ai_dot_aj) # tanimoto similarity
                heappush( min_heap, (1-S,edge_pair) )
    return [ heappop(min_heap) for i in range(len(min_heap)) ] # return ordered edge pairs









#if __name__ == '__main__':
#    # build option parser:
#    class MyParser(OptionParser):
#        def format_epilog(self, formatter):
#            return self.epilog
#    
#    usage = "usage: python %prog [options] filename"
#    description = """The link communities method of Ahn, Bagrow, and Lehmann, Nature, 2010:
#    www.nature.com/nature/journal/v466/n7307/full/nature09182.html (doi:10.1038/nature09182)
#    """
#    epilog = """
#    
#Input:
#  An edgelist file where each line represents an edge:
#    node_i <delimiter> node_j <newline>
#  if unweighted, or
#    node_i <delimiter> node_j <delimiter> weight_ij <newline>
#  if weighted.
#    
#Output: 
#  Three text files with extensions .edge2comm.txt, .comm2edges.txt,
#  and .comm2nodes.txt store the communities.
# 
#  edge2comm, an edge on each line followed by the community
#  id (cid) of the edge's link comm:
#    node_i <delimiter> node_j <delimiter> cid <newline>
#  
#  comm2edges, a list of edges representing one community per line:
#    cid <delimiter> ni,nj <delimiter> nx,ny [...] <newline>
#
#  comm2nodes, a list of nodes representing one community per line:
#    cid <delimiter> ni <delimiter> nj [...] <newline>
#  
#  The output filename contains the threshold at which the dendrogram
#  was cut, if applicable, or the threshold where the maximum
#  partition density was found, and the value of the partition 
#  density.
#  
#  If no threshold was given to cut the dendrogram, a file ending with
#  `_thr_D.txt' is generated, containing the partition density as a
#  function of clustering threshold.
#"""
#    parser = MyParser(usage, description=description,epilog=epilog)
#    parser.add_option("-d", "--delimiter", dest="delimiter", default="\t",
#                      help="delimiter of input & output files [default: tab]")
#    parser.add_option("-t", "--threshold", dest="threshold", type="float", default=None,
#                      help="threshold to cut the dendrogram (optional)")
#    parser.add_option("-w", "--weighted", dest="is_weighted", action="store_true", default=False,
#                    help="is the network weighted?")
#                    
#    # parse options:
#    (options, args) = parser.parse_args()
#    if len(args) != 1:
#        parser.error("incorrect number of arguments")
#    delimiter = options.delimiter
#    if delimiter == '\\t':
#        delimiter = '\t'
#    threshold   = options.threshold
#    is_weighted = options.is_weighted
#    
#    
#    print "# loading network from edgelist..."
#    basename = os.path.splitext(args[0])[0]
#    if is_weighted:
#        adj,edges,ij2wij = read_edgelist_weighted(args[0], delimiter=delimiter)
#    else: 
#        adj,edges        = read_edgelist_unweighted(args[0], delimiter=delimiter)
#    
#    
#    # run the method:
#    if threshold is not None:
#        if is_weighted:
#            edge2cid,D_thr = HLC( adj,edges ).single_linkage( threshold, w=ij2wij )
#        else:
#            edge2cid,D_thr = HLC( adj,edges ).single_linkage( threshold )
#        print "# D_thr = %f" % D_thr
#        write_edge2cid( edge2cid,"%s_thrS%f_thrD%f" % (basename,threshold,D_thr), delimiter=delimiter )
#    else:
#        if is_weighted:
#            edge2cid,S_max,D_max,list_D = HLC( adj,edges ).single_linkage( w=ij2wij )
#        else:
#            edge2cid,S_max,D_max,list_D = HLC( adj,edges ).single_linkage()
#        f = open("%s_thr_D.txt" % basename,'w')
#        for s,D in list_D:
#            print >>f, s, D
#        f.close()
#        print "# D_max = %f\n# S_max = %f" % (D_max,S_max)
#        write_edge2cid( edge2cid,"%s_maxS%f_maxD%f" % (basename,S_max,D_max), delimiter=delimiter )

