# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-01-21 13:00
# modified : 2013-04-09 08:38
"""
Tree-based detection of borrowings in lexicostatistical wordlists.
"""

__author_="Johann-Mattis List"
__date__="2013-04-09"


# basic imports
import os
import json
import itertools

# thirdparty imports
import numpy as np
import networkx as nx
import scipy.stats as sps
import numpy.linalg as linalg

# import error classes
from ...check.exceptions import *
from ...check.messages import *

# mpl is only used for specific plots, we can therefor make a safe import
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    ThirdPartyModuleError('matplotlib').warning()

# import the geoplot module
try:
    import mpl_toolkits.basemap as bmp
except ImportError:
    ThirdPartyModuleError('basemap').warning()

from .polygon import getConvexHull

# lingpy imports
from ...thirdparty import cogent as cg
from ...convert.gml import *
from ...basic import Wordlist
from ...read.csv import csv2dict,csv2list

class TreBor(Wordlist):
    """
    Basic class for calculations using the TreBor method.

        
    Parameters
    ----------
    dataset : string
        Name of the dataset that shall be analyzed.
    tree : {None, string}
        Name of the tree file.
    paps : string (default="pap")
        Name of the column that stores the specific cognate IDs consisting
        of an arbitrary integer key and a key for the concept.
    cognates : string (default="cogid")
        Name of the column that stores the general cognate ids.
    verbose : bool (default=False)
        Handle verbose output.
    tree_calc : {'neighbor','upgma'} (default='neighbor')
        Select the algorithm to be used for the tree calculation if no tree is
        passed with the file.

    """

    # XXX generally: find a way to check whether a dataset was already loaded,
    # XXX otherwise it takes too long a time to recalculate everything
    
    def __init__(
            self,
            dataset,
            tree = None,
            paps = 'pap',
            cognates = 'cogid',
            verbose = False,
            tree_calc = 'neighbor',
            **keywords
            ):
        # TODO check for keywords, allow to load trees, etc.

        # store the name of the dataset and the identifier for paps
        if dataset.endswith('.csv'):
            self.dataset = dataset.replace('.csv','')
        else:
            self.dataset = dataset
        self._pap_string = paps

        # open csv-file of the data and store it as a word list attribute
        Wordlist.__init__(self,self.dataset+'.csv',row='concept',col='doculect')
        
        #self.Wordlist = Wordlist(dataset+'.csv')
        self = self

        if verbose: print("[i] Loaded the wordlist file.")

        # check for glossid
        if 'glid' not in self.entries:
            self._gl2id = dict(
                    zip(
                        self.rows,
                        [i+1 for i in range(len(self.rows))]
                        )
                    )
            self._id2gl = dict([(b,a) for a,b in self._gl2id.items()])

            f = lambda x: self._gl2id[x]

            self.add_entries(
                    'glid',
                    'concept',
                    f
                    )

        # check for paps as attribute in the wordlist
        if paps not in self.entries:
            
            # define the function for conversion
            f = lambda x,y: "{0}:{1}".format(x[y[0]],x[y[1]])
            self.add_entries(
                    paps,
                    cognates+',glid',
                    f
                    )

            if verbose: print("[i] Created entry PAP.")
        
        # get the paps and the etymological dictionary
        if not hasattr(self,'paps'):
            self.paps = self.get_paps(ref=paps,missing=-1)
            self.etd = self.get_etymdict(ref=paps)

        if verbose: print("[i] Created the PAP matrix.")

        # get a list of concepts corresponding to the cogs and get the
        # singletons to be excluded from the calculation
        if not hasattr(self,'singletons'):
            tmp = self.get_etymdict(ref=paps,entry='concept')
            
            # a dictionary with pap-key as key and concept as value
            self.pap2con = {}

            # list stores the singletons
            self.singletons = []

            for key in self.paps:
                
                # get the names of the concepts
                concept_list = [k for k in tmp[key] if k != 0]
                concept = concept_list[0][0]
                self.pap2con[key] = concept

                # check for singletons
                if sum([1 for p in self.paps[key] if p >= 1]) == 1:
                    self.singletons.append(key)
            
            # create a list of keys for faster access when iterating
            self.cogs = [k for k in self.pap2con if k not in self.singletons]

            if verbose: print("[i] Excluded singletons.")

        # Load the tree, if it is not defined, assume that the treefile has the
        # same name as the dataset
        if not tree and not hasattr(self,'tree'):
            # try to load the tree first
            try:
                self.tree = cg.LoadTree(dataset+'.tre')
            except:
                # create it otherwise
                self.calculate(
                        'tree',
                        cognates=cognates,
                        tree_calc=tree_calc,
                        verbose=verbose
                        )
                if verbose: print("[i] Tree-file was not found, creating it now...")
            # XXX TODO
        
        # if it is explicitly defined, try to load that file
        elif not hasattr(self,'tree'):
            self.tree = cg.LoadTree(tree)
            if verbose: print("[i] Loaded the tree.")
        else:
            pass



        ## get the taxa
        ##self.taxa = self.cols

        #if verbose: print("[i] Assigned the taxa.")
        
        # create a couple of further attributes
        for a in ['stats','gls','dists','graph']:
            if not hasattr(self,a):
                setattr(self,a,{})
    
    def _get_GLS_top_down(
            self,
            pap,
            mode = 1,
            verbose = False
            ):
        """
        Infer gain-loss scenario using the method by Dagan & Martin (2007).

        """
        # check for mode
        try:
            mode = int(mode)
        except:
            raise ValueError("[i] Mode should be an integer.")
        
        # get the list of nodes that are not missing
        taxa,paps = [],[]
        for i,taxon in enumerate(self.taxa):
            if pap[i] != -1:
                taxa += [taxon]
                paps += [pap[i]]

        # get the subtree of all taxa
        tree = self.tree.getSubTree(taxa)

        # get list of taxa where pap is 1
        presents = [self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1]

        # get the subtree containing all taxa that have positive paps
        tree = tree.lowestCommonAncestor(
                [
                    self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1
                    ]
                )

        if verbose: print("[i] Subtree is {0}.".format(str(tree)))

        # assign the basic (starting) values to the dictionary
        nodes = [t.Name for t in tree.tips()]

        if verbose: print("[i] Nodes are {0}.".format(','.join(nodes)))
        
        if mode == 1:
            return [(tree.Name,1)]
        
        # store the scenario
        scenario = []

        # make the queue
        queue = [[tree,1]]
        while queue:
            
            # get tree and counter from queue
            tmp_tree,counter = queue.pop(0)

            # break if counter exceeds the mode
            if counter >= mode:
                t = tmp_tree.lowestCommonAncestor([p for p in presents if p in
                    tmp_tree.getTipNames()])
                scenario += [(t.Name,1)]
            else:
                # get tip names for checking
                tmp_names = tmp_tree.getTipNames()
                
                # get the children
                childA,childB = tmp_tree.Children

                # get lowest common ancestor
                subA = childA.lowestCommonAncestor([p for p in presents if p in
                    childA.getTipNames()])
                subB = childB.lowestCommonAncestor([p for p in presents if p in
                    childB.getTipNames()])

                # check for tip names in subtrees
                if subA.Children:
                    subAnames = subA.getTipNames()
                else:
                    if childA.Name in presents:
                        subAnames = [childA.Name]
                    else:
                        subAnames = []
                if subB.Children:
                    subBnames = subB.getTipNames()
                else:
                    if childB.Name in presents:
                        subBnames = [childB.Name]
                    else:
                        subBnames = []

                commons = subAnames + subBnames

                # check for identity, if the tips are identical, stop the iteration
                if set(tmp_names) == set(commons):
                    scenario += [(tmp_tree.lowestCommonAncestor(presents).Name,1)]

                # if they are not, append the subtrees to the queue
                else:
                    if subA.Children:
                        queue += [[subA,counter+1]]
                    else:
                        scenario += [(subA.Name,1)]
                    
                    if subB.Children:
                        queue += [[subB,counter+1]]
                    else:
                        scenario += [(subB.Name,1)]

        # TODO fill the scenario with gaps
        output = []
        d = {}
        for i,taxon in enumerate(taxa):
            if paps[i] >= 1:
                d[taxon] = 1
            else:
                d[taxon] = 0

        for s in scenario:
            
            # append scenario to output
            output += [s]

            # get the subtree
            subtree = tree.getNodeMatchingName(s[0])

            # check whether subtree is leave
            if not subtree.Children:
                pass
            else:
                # order the internal nodes according to the number of their leaves
                ordered_nodes = sorted(
                        subtree.nontips()+[subtree],key=lambda x:len(x.tips())
                        )

                # start bottom-up
                for node in ordered_nodes:
                    nA,nB = node.Children

                    # get the states
                    stateA = d[nA.Name]
                    stateB = d[nB.Name]

                    # compare the states
                    if stateA == 1 and stateB == 1:
                        d[node.Name] = 1
                    elif stateA == 0 and stateB == 0:
                        d[node.Name] = 0
                    else:
                        d[node.Name] = 1
                        if stateA == 0:
                            output += [(nA.Name,0)]
                        elif stateB == 0:
                            output += [(nB.Name,0)]

        return output

        
    def _get_GLS(
            self,
            pap,
            mode = 'w',
            r = (1,1),
            verbose = False
            ):
        """
        Calculate a gain-loss scenario (GLS) for a given PAP.

        """
        # make a dictionary that stores the scenario
        d = {}
        
        # get the list of nodes that are not missing
        taxa,paps = [],[]
        for i,taxon in enumerate(self.taxa):
            if pap[i] != -1:
                taxa += [taxon]
                paps += [pap[i]]

        # get the subtree of all taxa
        tree = self.tree.getSubTree(taxa)

        # get the subtree containing all taxa that have positive paps
        tree = tree.lowestCommonAncestor(
                [
                    self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1
                    ]
                )

        if verbose: print("[i] Subtree is {0}.".format(str(tree)))

        # assign the basic (starting) values to the dictionary
        nodes = [t.Name for t in tree.tips()]

        if verbose: print("[i] Nodes are {0}.".format(','.join(nodes)))

        # get the first state of all nodes and store the state in the
        # dictionary. note that we start from two distinct scenarios: one
        # assuming single origin at the root where all present states in the
        # leave are treated as retentions, and one assuming multiple origins,
        # where all present states in the leaves are treated as origins
        for node in nodes:
            idx = taxa.index(node)
            if paps[idx] >= 1:
                state = 1
            else:
                state = 0
            d[node] = [(state,[])]

        # return simple scenario, if the group is single-origin
        if sum([d[node][0][0] for node in nodes]) == len(nodes):
            return [(tree.Name,1)]

        # order the internal nodes according to the number of their leaves
        ordered_nodes = sorted(
                tree.nontips()+[tree],key=lambda x:len(x.tips())
                )

        # calculate the general restriction value (maximal weight). This is roughly
        # spoken simply the minimal value of either all events being counted as
        # origins (case 1) or assuming origin of the character at the root and
        # counting all leaves that lost the character as single loss events (case
        # 2). In case two, the first gain of the character has to be added
        # additionally
        if mode == 'w':
            RST = min(paps.count(1) * r[0], paps.count(0) * r[1] + r[0])
        elif mode == 'r':
            RST = r

        # join the nodes successively
        for i,node in enumerate(ordered_nodes):
            if verbose: print(node.Name)
            
            # when dealing with multifurcating trees, we have to store all
            # possible scenarios, i.e. we need to store the crossproduct of all
            # scenarios

            # get the names of the children of the nodes
            names = [x.Name for x in node.Children]

            # get the nodes with their states from the dictionary
            tmp_nodes = [d[x.Name] for x in node.Children]

            # get the cross-product of the stuff
            crossp = itertools.product(*tmp_nodes)
            
            newNodes = []
            
            # combine the histories of the items if all have the same value,
            # therefore, we first get the states in a simple list
            for cross in crossp:
                states = [x[0] for x in cross]
                stories = [x[1] for x in cross]
                states_sum = sum(states)
                states_len = len(states)
                
                # combine the histories
                new_stories = []
                for x in stories:
                    new_stories += x
                
                # if states are identical and point to gain / presence of chars
                if states_sum == states_len: 
                    # add the histories to the queue only if their weight is
                    # less or equal to the maxWeight
                    gl = [k[1] for k in new_stories]+[1]
                    
                    if mode == 'w':
                        weight = gl.count(1) * r[0] + gl.count(0) * r[1]
                    else:
                        weight = gl.count(1) + 1 # we need to add 1 here

                    if weight <= RST:
                        newNodes.append((1,new_stories))
                
                # if states are identical and point to absence of chars
                elif states_sum == 0:
                    gl = [k[1] for k in new_stories]
                    
                    if mode == 'w':
                        weight = gl.count(1) * r[0] + gl.count(0) * r[1]
                    else:
                        weight = gl.count(1)
                        
                    if weight <= RST:
                        newNodes.append((0,new_stories))

                # if the states are not identical, we check for both scenarios
                else:
                    # first scenario (tmpA) assumes origin, that is, for each node
                    # that has a 1, we add an origin to new_stories, same is
                    # for loss scenario (tmpB)
                    tmpA = [x for x in new_stories]
                    tmpB = [x for x in new_stories]
                    for c,state in enumerate(states):
                        if state == 1:
                            tmpA += [(names[c],1)]
                        else:
                            tmpB += [(names[c],0)]

                    # get the vectors to make it easier to retrieve the number
                    # of losses and gains
                    glA = [k[1] for k in tmpA]
                    glB = [k[1] for k in tmpB] + [1] # don't forget adding 1 origin

                    # check the gain-loss scores
                    if mode == 'w':
                        weightA = glA.count(1) * r[0] + glA.count(0) * r[1]
                        weightB = glB.count(1) * r[0] + glB.count(0) * r[1]
                    else:
                        weightA = glA.count(1)
                        weightB = glB.count(1)

                    newNodeA = (0,tmpA)
                    newNodeB = (1,tmpB)
                    
                    # check for additional gains in the gain-scenario,
                    # according to the current model, we don't allow for one
                    # character to be gained twice along a branch, i.e. by an
                    # ancestor, then get lost, and than be gained anew
                    if 1 in [k[1] for k in tmpB]:
                        noB = True
                    else:
                        noB = False

                    if weightA <= RST:
                        newNodes += [newNodeA]
                    if weightB <= RST and not noB:
                        newNodes += [newNodeB]
                        
                d[node.Name] = newNodes
                if verbose: print("nodelen",len(d[node.Name]))
        
        # try to find the best scenario by counting the ratio of gains and losses.
        # the key idea here is to reduce the number of possible scenarios according
        # to a given criterion. We choose the criterion of minimal changes as a
        # first criterion to reduce the possibilities, i.e. we weight both gains
        # and losses by 1 and select only those scenarios where gains and losses
        # sum up to a minimal number of gains and losses. This pre-selection of
        # scenarios can be further reduced by weighting gains and losses
        # differently. So in a second stage we choose only those scenarios where
        # there is a minimal amount of gains. 
        
        if verbose: print(len(d[tree.Name]))

        # convert the specific format of the d[tree.Name] to simple format
        gls_list = []
        for first,last in d[tree.Name]:
            if first == 1:
                gls_list.append([(tree.Name,first)]+last)
            else:
                gls_list.append(last)

        # the tracer stores all scores
        tracer = []

        for i,line in enumerate(gls_list):
            
            # calculate gains and losses
            gains = sum([1 for x in line if x[1] == 1])
            losses = sum([1 for x in line if x[1] == 0])

            # calculate the score
            if mode == 'w':
                score = r[0] * gains + r[1] * losses
            else:
                score = gains + losses

            # append it to the tracer
            tracer.append(score)
        
        # get the minimum score
        minScore = min(tracer)

        if mode == 'w':

            # return the minimal indices, sort them according to the number of
            # gains inferred, thereby pushing gains to the root, similar to
            # Mirkin's (2003) suggestion
            return sorted(
                    [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore],
                    key = lambda x:sum([i[1] for i in x])
                    )[0]
        
        # push gains down to the root as suggested by Mirkin 2003
        minimal_gains = [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore]
        
        best_scenario = None
        old_length_of_tips = len(self.taxa) + 1

        for i,line in enumerate(minimal_gains):
            
            # calculate number of tips for the gains of a given scenario
            new_length_of_tips = 0
            for taxon,state in line:
                if state == 1:
                    new_length_of_tips += len(
                            self.tree.getNodeMatchingName(taxon).getTipNames()
                            )
            if new_length_of_tips < old_length_of_tips:
                old_length_of_tips = new_length_of_tips
                best_scenario = i
        return minimal_gains[best_scenario]

    def get_GLS(
            self,
            mode = 'weighted',
            ratio = (1,1),
            restriction = 3,
            output_gml = False,
            output_plot = False,
            verbose = False,
            tar = False
            ):
        """
        Create gain-loss-scenarios for all non-singleton paps in the data.

        Parameters
        ----------
        mode : string (default="weighted")
            Select between "weighted" and "restriction".
        ratio : tuple (default=(1,1))
            If "weighted" mode is selected, define the ratio between the
            weights for gains and losses.
        restriction : int (default=3)
            If "restriction" is selected as mode, define the maximal number of
            gains.
        output_gml : bool (default=False)
            If set to c{True}, the decisions for each GLS are stored in a
            separate file in GML-format.
        tar : bool (default=False)
            If set to c{True}, the GML-files will be added to a compressed tar-file.

        """
        if mode not in ['weighted','w','r','restriction','t','topdown']:
            raise ValueError("[!] The mode {0} is not available".format(mode))

        # define alias for mode
        if mode in ['w','weighted']:
            mode = 'weighted'
        elif mode in ['r','restriction']:
            mode = 'restriction'
        else:
            mode = 'topdown'

        # create a named string for the mode
        if mode == 'weighted':
            glm = 'w-{0[0]}-{0[1]}'.format(ratio)
        elif mode == 'restriction':
            glm = 'r-{0}'.format(restriction)
        elif mode == 'topdown':
            glm = 't-{0}'.format(restriction)
        
        # create statistics for this run
        self.stats[glm] = {}

        # store the statistics
        self.stats[glm]['mode'] = mode
        self.stats[glm]['dataset'] = self.dataset

        # attribute stores all gls for each cog
        self.gls[glm] = {}

        for cog in self.cogs:
            if verbose: print("[i] Calculating GLS for COG {0}...".format(cog),end="")
            if mode == 'weighted':
                gls = self._get_GLS(
                        self.paps[cog],
                        r = ratio,
                        mode = 'w'
                        )

            if mode == 'restriction':
                gls = self._get_GLS(
                        self.paps[cog],
                        r = restriction,
                        mode = 'r'
                        )

            if mode == 'topdown':
                gls = self._get_GLS_top_down(
                        self.paps[cog],
                        mode = restriction
                        )

            noo = sum([t[1] for t in gls])
            
            self.gls[glm][cog] = (gls,noo)


            # attend scenario to gls
            if verbose: print(" done.")
        if verbose: print("[i] Successfully calculated Gain-Loss-Scenarios.")
 
        # write the results to file
        # make the folder for the data to store the stats
        folder = self.dataset+'_trebor'
        try:
            os.mkdir(folder)
        except:
            pass       
        
        # if output of gls is chosen, load the gml-graph
        if output_gml:

            # make the directory for the files
            try:
                os.mkdir(folder+'/gml')
            except:
                pass

            # make next directory
            try:
                os.mkdir(
                        folder+'/gml/'+'{0}-{1}'.format(
                            self.dataset,
                            glm
                            )
                        )
            except:
                pass

            # make the folder for png
            try:
                os.mkdir(
                        folder+'/gml/'+'{0}-{1}-figures'.format(
                            self.dataset,
                            glm
                            )
                        )
            except:
                pass

            # load the graph
            self.gml = nx.read_gml(self.dataset+'.gml')

            # store the graph
            for cog in self.cogs:
                gls = self.gls[glm][cog][0]
                g = gls2gml(
                        gls,
                        self.gml,
                        self.tree,
                        filename = folder+'/gml/{0}-{1}/{2}'.format(
                            self.dataset,
                            glm,
                            cog
                            ),
                        )

                # if plot of gml is chose
                if output_plot:
                    nodes = []
                    
                    for n,d in g.nodes(data=True):
                        x = d['graphics']['x']
                        y = d['graphics']['y']
                        f = d['graphics']['fill']
                        o = d['origin']
                        l = d['label']
                        
                        nodes.append((x,y,f,o,l))

                    edges = []
                    for a,b,d in g.edges(data=True):
                    
                        xA = g.node[a]['graphics']['x']
                        xB = g.node[b]['graphics']['x']
                        yA = g.node[a]['graphics']['y']
                        yB = g.node[b]['graphics']['y']
                    
                        edges += [(xA,xB,yA,yB)]
                    
                    #mpl.rc('text',usetex=keywords['usetex'])
                    fig = plt.figure()
                    figsp = fig.add_subplot(111)
                    ax = plt.axes(frameon=False)
                    plt.xticks([])
                    plt.yticks([])
                    
                    plt.axis('equal')
                    
                    for xA,xB,yA,yB in edges:
                    
                        plt.plot(
                                [xA,xB],
                                [yA,yB],
                                '-',
                                color='black',
                                linewidth=5
                                )
                        plt.plot(
                                [xA,xB],
                                [yA,yB],
                                '-',
                                color='0.2',
                                linewidth=4
                                )
                    for x,y,f,o,l in nodes:
                        if f == '#000000':
                            c = '#ffffff'
                        else:
                            c = '#000000'
                        if o == 1:
                            size = 20
                        else:
                            size = 10
                        if l.startswith('edge') or l.startswith('root'):
                            plt.plot(x,y,'o',markersize=size,color=f)
                        else:
                            plt.text(
                                    x,
                                    y,
                                    l,
                                    horizontalalignment='center',
                                    verticalalignment='center',
                                    size=8,fontweight='bold',color=c,backgroundcolor=f)
                    
                    #plt.subplots_adjust(left=0.02,right=0.98,top=0.98,bottom=0.02)
                    plt.savefig(folder+'/gml/{0}-{1}-figures/{2}.png'.format(
                        self.dataset,
                        glm,
                        cog
                        ))
                    plt.clf()

            # if tar is chosen, put it into a tarfile
            if tar:
                os.system(
                        'cd {0}_trebor/gml/ ; tar -pczf {0}-{1}.tar.gz {0}-{1}; cd ..; cd ..'.format(
                            self.dataset,
                            glm
                            )
                        )
                os.system('rm {0}_trebor/gml/{0}-{1}/*.gml'.format(self.dataset,glm))
                os.system('rmdir {0}_trebor/gml/{0}-{1}'.format(self.dataset,glm))


        # store some statistics as attributes
        self.stats[glm]['ano'] = sum(
                [v[1] for v in self.gls[glm].values()]
                ) / len(self.gls[glm])
        self.stats[glm]['mno'] = max([v[1] for v in self.gls[glm].values()])
        self.stats[glm]['ratio'] = ratio 
        self.stats[glm]['restriction'] = restriction

        # store statistics and gain-loss-scenarios in textfiles
        # create folder for gls-data
        try:
            os.mkdir(folder+'/gls')
        except:
            pass
        
        if verbose: print("[i] Writing GLS data to file... ",end="")
        
        # write gls-data to folder
        f = open(folder+'/gls/{0}-{1}.gls'.format(self.dataset,glm),'w')
        f.write('PAP\tGainLossScenario\tNumberOfOrigins\n')
        for cog in sorted(self.gls[glm]):
            gls,noo = self.gls[glm][cog]
            f.write(
                    "{0}\t".format(cog)+','.join(
                        ["{0}:{1}".format(a,b) for a,b in gls]
                        ) + '\t'+str(noo)+'\n'
                    )
        f.close()
        if verbose: print("done.")

        
        # print out average number of origins
        if verbose: print("[i] Average Number of Origins: {0:.2f}".format(self.stats[glm]['ano']))

        # write statistics to stats file
        try:
            os.mkdir(folder+'/stats')
        except:
            pass

        f = open(folder+'/stats/{0}-{1}'.format(self.dataset,glm),'w')
        f.write('Number of PAPs (total): {0}\n'.format(len(self.paps)))
        f.write('Number of PAPs (non-singletons): {0}\n'.format(len(self.gls[glm])))
        f.write('Number of Singletons: {0}\n'.format(len(self.singletons)))
        f.write('Average Number of Origins: {0:.2f}\n'.format(self.stats[glm]['ano']))
        f.write('Maximum Number of Origins: {0}\n'.format(self.stats[glm]['mno']))
        f.write('Mode: {0}\n'.format(mode))
        if mode == 'weighted':
            f.write('Ratio: {0[0]} / {0[1]}\n'.format(ratio))
        elif mode == 'restriction':
            f.write('Restriction: {0}\n'.format(restriction))

        f.close()

        return

    def get_CVSD(
            self,
            verbose = False
            ):
        """
        Calculate the Contemporary Vocabulary Size Distribution (CVSD).

        """
        # define taxa and concept as attribute for convenience
        taxa = self.taxa
        concepts = self.concept #XXX do we need this? XXX

        # calculate vocabulary size
        forms = []
        meanings = []
        for taxon in taxa:
            f = [x for x in set(
                self.get_list(col=taxon,entry=self._pap_string,flat=True)
                ) if x in self.cogs
                ]
            m = set([x.split(':')[1] for x in f])
            forms += [len(f)]
            meanings += [len(m)]
        
        # store the stuff as an attribute
        self.dists['contemporary'] = [x for x,y in zip(forms,meanings)] # XXX

        if verbose: print("[i] Calculated the distributions for contemporary taxa.")
        
        return 

    def get_AVSD(
            self,
            glm,
            verbose = False,
            write = False
            ):
        """
        Function retrieves all paps for ancestor languages in a given tree.
        """
        # define concepts for convenience
        concepts = self.concepts # XXX do we need this? XXX
        
        # get all internal nodes, i.e. the nontips and also the root
        nodes = ['root'] + sorted(
                [node.Name for node in self.tree.nontips()],
                key=lambda x: len(self.tree.getNodeMatchingName(x).tips()),
                reverse = True
                )

        # retrieve scenarios
        tmp = sorted([(a,b,c) for a,(b,c) in self.gls[glm].items()])
        cog_list = [t[0] for t in tmp]
        gls_list = [t[1] for t in tmp]
        noo_list = [t[2] for t in tmp]

        # create a list that stores the paps
        paps = [[0 for i in range(len(nodes))] for j in range(len(cog_list))]

        # iterate and assign values
        for i,cog in enumerate(cog_list):
            
            # sort the respective gls
            gls = sorted(
                    gls_list[i],
                    key = lambda x: len(self.tree.getNodeMatchingName(x[0]).tips()),
                    reverse = True
                    )

            # retrieve the state of the root
            if gls[0][1] == 1 and gls[0][0] == 'root':
                state = 1
            else:
                state = 0

            # assign the state of the root to all nodes
            paps[i] = [state for node in nodes]

            # iterate over the gls and assign the respective values to all
            # children
            for name,event in gls:
                if event == 1:
                    this_state = 1
                else:
                    this_state = 0

                # get the subtree nodes
                sub_tree_nodes = [node.Name for node in
                        self.tree.getNodeMatchingName(name).nontips()]

                # assign this state to all subtree nodes
                for node in sub_tree_nodes:
                    paps[i][nodes.index(node)] = this_state

        # get number of forms and number of meanings
        # extract cogs instead of numbers, XXX this can actually be done in the
        # step before, it's just for testing at the moment
        for i,cog in enumerate(cog_list):
            for j,t in enumerate(paps[i]):
                if t == 1:
                    paps[i][j] = cog
                else:
                    pass
        
        # get forms and meanings
        forms = []
        meanings = []
        for i in range(len(paps[0])):
            f = set([x[i] for x in paps if x[i] != 0])
            m = set([x[i].split(':')[1] for x in paps if x[i] != 0])
            forms += [len(f)]
            meanings += [len(m)]

        # store the number of forms as an attribute
        self.dists[glm] = [x for x,y in zip(forms,meanings)] # XXX

        # write results of analysis to file 
        # XXX adjust later
        if write:
            cogs = [k[self.header['pap']] for k in self._data.values()]
            protos = [k[self.header['proto']] for k in
                    self._data.values()]

            cogs2proto = dict(zip(cogs,protos))
            f = open('ancestral_paps.txt','w')
            for i,n in enumerate(nodes):
                for j,p in enumerate(paps):
                    #print(paps[i])
                    c = paps[j][i]
                    if c != 0:
                        m = self.pap2con[c]
                        p = cogs2proto[c]
                        
                        if n != "root":
                            node = self.tree.getNodeMatchingName(n)
                            node = str(node).replace('(','').replace(')','').replace(',','-')
                        else:
                            node = n
                        f.write('{0}\t{1}\t{2}\t{3}\n'.format(
                            node,
                            c,
                            m,
                            p
                            ))
            f.close()


        if verbose: print("[i] Calculated the distributions for ancestral taxa.")

        return

    def get_IVSD(
            self,
            verbose = False,
            output_gml = False,
            output_plot = False,
            tar = True
            ):
        """
        Calculate VSD on the basis of each item.

        """

        # define concepts and taxa for convenience
        concepts = self.concepts
        taxa = self.taxa

        # get all internal nodes, i.e. the nontips and also the root
        nodes = ['root'] + sorted(
                [node.Name for node in self.tree.nontips()],
                key=lambda x: len(self.tree.getNodeMatchingName(x).tips()),
                reverse = True
                )
        
        # make dictionary that stores the best models for each cognate set
        best_models = {}

        # make array for all nodes and a dict for the scenarios
        all_avsd = [0 for node in nodes]
        scenarios = {}

        # iterate over concepts
        for concept in concepts:

            # get paps
            tmp = self.get_dict(row=concept,entry=self._pap_string)

            # add to list if value is missing
            for taxon in taxa:
                if taxon not in tmp:
                    tmp[taxon] = []

            # calculate distribution for contemporary taxa
            cvsd = [len([i for i in tmp[j] if i in self.cogs]) for j in taxa]

            # calculate ancestral dists, get all paps first
            pap_set = [i for i in set(
                    self.get_list(
                        row=concept,
                        entry=self._pap_string,
                        flat=True
                        )
                    ) if i in self.cogs]

            # get the models
            models = sorted(list(self.gls.keys()))

            # get the scenarios
            avsd_list = []
            for idx,glm in enumerate(models):
                avsd_list += [[0 for node in nodes]]
                for pap in pap_set:
                    gls,noo = self.gls[glm][pap]

                    # sort the gls
                    gls = sorted(
                            gls,
                            key = lambda x: len(self.tree.getNodeMatchingName(x[0]).tips()),
                            reverse = True
                            )

                    # retrieve the state of the root
                    if gls[0][1] == 1 and gls[0][0] == 'root':
                        state = 1
                    else:
                        state = 0

                    # assign the state of the root to all nodes in tmp
                    tmp = [state for node in nodes]

                    # iterate over the gls and assign the respective values to all
                    # children
                    for name,event in gls:
                        if event == 1:
                            this_state = 1
                        else:
                            this_state = 0

                        # get the subtree nodes
                        sub_tree_nodes = [node.Name for node in
                                self.tree.getNodeMatchingName(name).nontips()]

                        # assign this state to all subtree nodes
                        for node in sub_tree_nodes:
                            tmp[nodes.index(node)] = this_state

                    # add the values to the avsd_list
                    avsd_list[-1] = [a+b for a,b in zip(avsd_list[-1],tmp)]
            
            # check for correctness
            if verbose: print(concept,avsd_list[-1],models[-1])
            
            # calculate best distribution
            zp_vsd = []
            cvsd_set = set(cvsd)
            for avsd in avsd_list:
                if len(cvsd_set) == 1 and set(avsd):
                    zp_vsd.append((0,1.0))
                else:
                    vsd = sps.mannwhitneyu(
                            cvsd,
                            avsd
                            )

                zp_vsd.append(vsd)
            
            # extract p-values
            p_vsd = [p for z,p in zp_vsd]
            maxP = max(p_vsd)
            maxIdx = p_vsd.index(maxP)
            best_model = models[maxIdx]

            for p in pap_set:
                gls,noo = self.gls[best_model][p]
                best_models[p] = (best_model,noo,maxP)
                scenarios[p] = (gls,noo)

            # add sum to general model
            all_avsd = [a+b for a,b in zip(avsd_list[maxIdx],all_avsd)]

        
        self.best_models = best_models
        #print(sum([n for m,n,o in best_models.values()]) / len(best_models))
        
        # append to distributions
        self.dists['mixed'] = all_avsd

        # append to available models
        self.gls['mixed'] = scenarios

        # write the results to file
        # make the folder for the data to store the stats
        folder = self.dataset+'_trebor'
        try:
            os.mkdir(folder)
        except:
            pass       
        
        # if output of gls is chosen, load the gml-graph
        if output_gml:

            # make the directory for the files
            try:
                os.mkdir(folder+'/gml')
            except:
                pass

            # make next directory
            try:
                os.mkdir(
                        folder+'/gml/'+'{0}-{1}'.format(
                            self.dataset,
                            "mixed"
                            )
                        )
            except:
                pass

            # make the folder for png
            try:
                os.mkdir(
                        folder+'/gml/'+'{0}-{1}-figures'.format(
                            self.dataset,
                            "mixed"
                            )
                        )
            except:
                pass

            # load the graph
            self.gml = nx.read_gml(self.dataset+'.gml')

            # store the graph
            for cog in self.cogs:
                gls = self.gls["mixed"][cog][0]
                g = gls2gml(
                        gls,
                        self.gml,
                        self.tree,
                        filename = folder+'/gml/{0}-{1}/{2}'.format(
                            self.dataset,
                            "mixed",
                            cog
                            ),
                        )

                # if plot of gml is chose
                if output_plot:
                    nodes = []
                    
                    for n,d in g.nodes(data=True):
                        x = d['graphics']['x']
                        y = d['graphics']['y']
                        f = d['graphics']['fill']
                        o = d['origin']
                        l = d['label']
                        
                        nodes.append((x,y,f,o,l))

                    edges = []
                    for a,b,d in g.edges(data=True):
                    
                        xA = g.node[a]['graphics']['x']
                        xB = g.node[b]['graphics']['x']
                        yA = g.node[a]['graphics']['y']
                        yB = g.node[b]['graphics']['y']
                    
                        edges += [(xA,xB,yA,yB)]
                    
                    #mpl.rc('text',usetex=keywords['usetex'])
                    fig = plt.figure()
                    figsp = fig.add_subplot(111)
                    ax = plt.axes(frameon=False)
                    plt.xticks([])
                    plt.yticks([])
                    
                    plt.axis('equal')
                    
                    for xA,xB,yA,yB in edges:
                    
                        plt.plot(
                                [xA,xB],
                                [yA,yB],
                                '-',
                                color='black',
                                linewidth=5
                                )
                        plt.plot(
                                [xA,xB],
                                [yA,yB],
                                '-',
                                color='0.2',
                                linewidth=4
                                )
                    for x,y,f,o,l in nodes:
                        if f == '#000000':
                            c = '#ffffff'
                        else:
                            c = '#000000'
                        if o == 1:
                            size = 20
                        else:
                            size = 10
                        if l.startswith('edge') or l.startswith('root'):
                            plt.plot(x,y,'o',markersize=size,color=f)
                        else:
                            plt.text(
                                    x,
                                    y,
                                    l,
                                    horizontalalignment='center',
                                    verticalalignment='center',
                                    size=8,fontweight='bold',color=c,backgroundcolor=f)
                    
                    #plt.subplots_adjust(left=0.02,right=0.98,top=0.98,bottom=0.02)
                    plt.savefig(folder+'/gml/{0}-{1}-figures/{2}.png'.format(
                        self.dataset,
                        "mixed",
                        cog
                        ))
                    plt.clf()

            # if tar is chosen, put it into a tarfile
            if tar:
                os.system(
                        'cd {0}_trebor/gml/ ; tar -pczf {0}-{1}.tar.gz {0}-{1}; cd ..; cd ..'.format(
                            self.dataset,
                            "mixed"
                            )
                        )
                os.system('rm {0}_trebor/gml/{0}-{1}/*.gml'.format(self.dataset,"mixed"))
                os.system('rmdir {0}_trebor/gml/{0}-{1}'.format(self.dataset,"mixed"))



        # store some statistics as attributes
        self.stats['mixed'] = {}
        self.stats['mode'] = 'mixed'
        self.stats['dataset'] = self.dataset
        self.stats['mixed']['ano'] = sum(
                [v[1] for v in self.gls['mixed'].values()]
                ) / len(self.gls['mixed'])
        self.stats['mixed']['mno'] = max([v[1] for v in self.gls['mixed'].values()])

        # store statistics and gain-loss-scenarios in textfiles
        # create folder for gls-data
        try:
            os.mkdir(folder+'/gls')
        except:
            pass
        
        if verbose: print("[i] Writing GLS data to file... ",end="")
        
        # write gls-data to folder
        f = open(folder+'/gls/{0}-{1}.gls'.format(self.dataset,"mixed"),'w')
        f.write('PAP\tGainLossScenario\tNumberOfOrigins\n')
        for cog in sorted(self.gls["mixed"]):
            gls,noo = self.gls["mixed"][cog]
            f.write(
                    "{0}\t".format(cog)+','.join(
                        ["{0}:{1}".format(a,b) for a,b in gls]
                        ) + '\t'+str(noo)+'\n'
                    )
        f.close()
        if verbose: print("done.")

        return 

    def get_MLN(
            self,
            glm,
            threshold = 1,
            verbose = False,
            #colormap = None,
            degree = 100
            ):
        """
        Compute an Minimal Lateral Network for a given model.

        Parameters
        ----------
        glm : str 
            The dictionary key for the gain-loss-model.
        threshold : int (default=1)
            The threshold used to exclude edges.
        verbose : bool (default=False)
            Set to c{True} for verbose output.
        colormap : {None matplotlib.cm}
            A :py:class:`matplotlib.colormap` instance. If set to c{None}, this
            defaults to :py:class:`~matplotlib.cm.jet`.
        degree : int (default=100)
            The degree which is chosen for the projection of the tree layout.
        """

        #if not colormap:
        #    colormap = mpl.cm.jet
        
        # create the primary graph
        gPrm = nx.Graph()

        # make alias for tree and taxa for convenience
        taxa = self.taxa
        tree = self.tree

        # create the template graph XXX add fallback procedure
        try:
            gTpl = nx.read_gml(self.dataset+'.gml')
        except:
            
            # if no good topology is given, create it automatically, using
            # the radial layout function
            gTpl = radial_layout(str(self.tree),filename='',degree=degree)
            
            if verbose: print("[i] Calculated radial layout for the tree. ")
            
        # make alias for the current gls for convenience
        scenarios = self.gls[glm]

        # create dictionary for inferred lateral events
        ile = {}

        # create mst graph
        gMST = nx.Graph()

        # create out graph
        gOut = nx.Graph()

        # load data for nodes into new graph
        for node,data in gTpl.nodes(data=True):
            if data['label'] in taxa:
                data['graphics']['fill'] = '#ff0000'
                data['graphics']['type'] = 'rectangle'
                data['graphics']['w'] = 80.0
                data['graphics']['h'] = 20.0
            else:
                data['graphics']['type'] = 'ellipse'
                data['graphics']['w'] = 30.0
                data['graphics']['h'] = 30.0
                data['graphics']['fill'] = '#ff0000'
            gPrm.add_node(data['label'],**data)

        # load edge data into new graph
        for nodeA,nodeB,data in gTpl.edges(data=True):
            if 'graphics' not in data:
                data['graphics'] = {}
            data['graphics']['width'] = 10.0
            data['graphics']['fill'] = '#000000'
            data['label'] = 'vertical'

            gPrm.add_edge(
                    gTpl.node[nodeA]['label'],
                    gTpl.node[nodeB]['label'],
                    **data
                    )

        # start to assign the edge weights
        for cog,(gls,noo) in scenarios.items():
            
            # get the origins
            oris = [x[0] for x in gls if x[1] == 1]

            # connect origins by edges
            for i,oriA in enumerate(oris):
                for j,oriB in enumerate(oris):
                    if i < j:
                        try:
                            gPrm.edge[oriA][oriB]['weight'] += 1
                        except:
                            gPrm.add_edge(
                                    oriA,
                                    oriB,
                                    weight=1
                                    )

        # verbose output
        if verbose: print("[i] Calculated primary graph.")
        
        # verbose output
        if verbose: print("[i] Inferring lateral edges...")
            
        # create MST graph
        gMST = nx.Graph()

        for cog,(gls,noo) in scenarios.items():
            
            ile[cog] = []

            # get the origins
            oris = [x[0] for x in gls if x[1] == 1]

            # create a graph of weights
            gWeights = nx.Graph()
            
            # iterate over nodes
            for i,nodeA in enumerate(oris):
                for j,nodeB in enumerate(oris):
                    if i < j:
                        w = gPrm.edge[nodeA][nodeB]['weight']
                        gWeights.add_edge(
                                nodeA,
                                nodeB,
                                weight=w
                                )
            
            # if the graph is not empty
            if gWeights:

                # check for identical weights and change them according to
                # tree-distance
                tmp_weights = {}
                for a,b,d in gWeights.edges(data=True):
                    try:
                        tmp_weights[int(d['weight'])] += [(a,b)]
                    except:
                        tmp_weights[int(d['weight'])] = [(a,b)]
                
                # check for identical weights and calculate the tree distance
                for w in tmp_weights:
                    elist = tmp_weights[w]
                    
                    # check whether there are more identical weights
                    if len(elist) > 1:
                        
                        # if so, order all stuff according to branch length
                        branches = []
                        for a,b in elist:
                            branch_distance = len(self.tree.getConnectingEdges(a,b))
                            branches += [(a,b,branch_distance)]

                        # now change the weights according to the order
                        scaler = 1 / len(branches)
                        minus = 1 - scaler
                        branches = sorted(branches,key=lambda x:(x[2],x[1],x[0]),reverse=True)
                        for a,b,d in branches:
                            gWeights.edge[a][b]['weight'] += minus
                            minus -= scaler
                
                # change maximum weights to distance weights
                for a,b,d in sorted(gWeights.edges(data=True),key=lambda x:x[2]['weight']):
                    w = d['weight']
                    gWeights.edge[a][b]['weight'] = int(1000 / w)
                
                # calculate the MST
                mst = nx.minimum_spanning_tree(gWeights,weight='weight')
                
                # assign the MST-weights to gMST
                for nodeA,nodeB in mst.edges():
                    try:
                        gMST.edge[nodeA][nodeB]['weight'] += 1
                        gMST.edge[nodeA][nodeB]['cogs'] += [cog]
                    except:
                        gMST.add_edge(
                                nodeA,
                                nodeB,
                                weight=1,
                                cogs=[cog]
                                )
                    ile[cog]+= [(nodeA,nodeB)]

        # load data for nodes into new graph
        for node,data in gTpl.nodes(data=True):
            if data['label'] in taxa:
                data['graphics']['fill'] = '#ff0000'
                data['graphics']['type'] = 'rectangle'
                data['graphics']['w'] = 80.0
                data['graphics']['h'] = 20.0
            else:
                data['graphics']['type'] = 'ellipse'
                data['graphics']['w'] = 30.0
                data['graphics']['h'] = 30.0
                data['graphics']['fill'] = '#ff0000'
            
            gOut.add_node(data['label'],**data)

        # load edge data into new graph
        for nodeA,nodeB,data in gTpl.edges(data=True):
            data['graphics']['width'] = 10.0
            data['graphics']['fill'] = '#000000'
            data['label'] = 'vertical'
            try:
                del data['graphics']['Line']
            except:
                pass
            gOut.add_edge(
                    gTpl.node[nodeA]['label'],
                    gTpl.node[nodeB]['label'],
                    **data
                    )
        
        # assign new edge weights
        for nodeA,nodeB,data in gMST.edges(data=True):
            w = data['weight']

            # get the color for the weight
            #color = mpl.colors.rgb2hex(colormap(cfunc[weights.index(w)]))

            data['graphics'] = {}
            #data['graphics']['fill'] = color
            #data['graphics']['width'] = w * scale
            data['cogs'] = ','.join([str(i) for i in data['cogs']])
            data['label'] = 'horizontal'

            # check for threshold
            if w >= threshold:
                try:
                    gOut.edge[nodeA][nodeB]
                except:
                    # add the data to the out-graph
                    gOut.add_edge(
                        nodeA,
                        nodeB,
                        **data
                        )
        # transfer node data

        # verbose output
        if verbose: print("[i] Writing graph to file...")

        # write the graph to file
        f = open(self.dataset+'_trebor/mln-'+glm+'.gml','w')
        for line in nx.generate_gml(gOut):
            f.write(line+'\n')
        f.close()

        # write the inferred borrowing events (ILS, inferred lateral event) 
        # between all taxa to file
        # verbose output
        if verbose: print("[i] Writing Inferred Lateral Events to file...")

        f = open(self.dataset+'_trebor/ile-'+glm+'.csv','w')
        for cog,events in ile.items():
            if events:
                f.write(
                        cog+'\t'+','.join(
                            ['{0}:{1}'.format(a,b) for a,b in events]
                            )+'\n'
                        )
        f.close()

        # create file name for node labels (cytoscape output)
        f = open(self.dataset+'_trebor/node.label.NA','w')
        f.write("node.label (class=java.lang.String)\n")
        for taxon in taxa:
            f.write('{0} = {1}\n'.format(taxon,taxon))
        f.close()

        # add gOut to graphattributes
        self.graph[glm] = gOut

        # write stats to file
        f = open(self.dataset+'_trebor/taxa-'+glm+'.stats','w')
        
        # get the degree
        nodes = tree.getNodeNames()

        dgr,wdgr = [],[]
        for taxon in nodes:
            
            horizontals = [g for g in gOut[taxon] if 'weight' in gOut[taxon][g]]
            
            dgr.append(len(horizontals))
            wdgr.append(sum([gOut[taxon][g]['weight'] for g in horizontals]))

        sorted_nodes = sorted(
                zip(nodes,dgr,wdgr),
                key=lambda x:x[1],
                reverse=True
                )
        for n,d,w in sorted_nodes:
            f.write(
                    '{0}\t{1}\t{2}\t{3}\n'.format(
                        n,
                        str(tree.getNodeMatchingName(n)),
                        d,
                        w
                        )
                    )
        f.close()

        if verbose: print("[i] Wrote node degree distributions to file.")

        # write edge distributions
        f = open(self.dataset+'_trebor/edge-'+glm+'.stats','w')
        edges = []
        edges = [g for g in gOut.edges(data=True) if 'weight' in g[2]]

        for nA,nB,d in sorted(
                edges,
                key=lambda x: x[2]['weight'],
                reverse = True
                ):
            f.write(
                    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                        nA,
                        nB,
                        d['weight'],
                        d['cogs'],
                        tree.getNodeMatchingName(nA),
                        tree.getNodeMatchingName(nB)
                        )
                    )
        f.close()
        if verbose: print("[i] Wrote edge-weight distributions to file.")
        
        # write specific links of taxa to file
        try:
            os.mkdir(self.dataset+'_trebor/taxa-'+glm)
        except:
            pass

        for taxon in self.taxa:
            f = open(self.dataset+'_trebor/taxa-'+glm+'/'+taxon+'.csv','w')
            keys = [n for n in gOut[taxon] if gOut[taxon][n]['label'] == 'horizontal']
            for key in sorted(keys,key=lambda x:gOut[taxon][x]['weight']):
                for cog in gOut[taxon][key]['cogs'].split(','):
                    tmp = [x for x in self.etd[cog] if x != 0]
                    idx = [x[0] for x in tmp][0]
                    concept = self[idx,'concept']

                    if 'proto' in self.entries:
                        proto = cog + '\t'+self[idx,'proto']
                    else:
                        proto = cog
                    if 'note' in self.entries:
                        proto += '\t'+self[idx,'note']

                    f.write('{0}\t{1}\t{2}\n'.format(
                        key,
                        proto,
                        concept
                        ))
            f.close()
        if verbose: pirnt("[i] Wrote list of edges per taxa to file.")

        return 

    def get_PDC(
            self,
            glm,
            verbose = False,
            **keywords
            ):
        """
        Calculate Patchily Distributed Cognates.
        """
        
        patchy = {}
        paps = []

        for key,(gls,noo) in self.gls[glm].items():

            # get the origins
            oris = [x[0] for x in gls if x[1] == 1]
            
            # get the tip-taxa for each origin
            tips = []

            # get the losses 
            tmp_loss = [x[0] for x in gls if x[1] == 0]
            losses = []
            for l in tmp_loss:
                losses += self.tree.getNodeMatchingName(l).getTipNames()

            for i,ori in enumerate(oris):
                tips += [
                        (
                            i+1,
                            [t for t in self.tree.getNodeMatchingName(
                                ori
                                ).getTipNames() if t not in losses]
                            )
                        ]

            # now, all set of origins with their tips are there, we store them
            # in the patchy dictionary, where each taxon is assigned the
            # numerical value of the given patchy dist
            patchy[key] = {}
            if len(tips) > 1:
                for i,tip in tips:
                    for taxon in tip:
                        patchy[key][taxon] = i
            else:
                for i,tip in tips:
                    for taxon in tip:
                        patchy[key][taxon] = 0
            
            paps.append((key,noo))
        
        if verbose: print("[i] Retrieved patchy distributions.")

        # get the index for the paps in the wordlist
        papIdx = self.header['pap']
        taxIdx = self._colIdx

        # create a dictionary as updater for the wordlist
        updater = {}
        for key in self:

            # get the taxon first
            taxon = self[key][taxIdx]

            # get the pap
            pap = self[key][papIdx]
            
            try:
                updater[key] = '{0}:{1}'.format(pap,patchy[pap][taxon])
            except KeyError:
                updater[key] = '{0}:{1}'.format(pap,0)

        # update the wordlist
        self.add_entries(
                'patchy',
                updater,
                lambda x:x
                )

        # write data to file
        # self.output('csv',filename=self.dataset+'_trebor/wl-'+glm)
        # XXX change later

        if verbose: print("[i] Updated the wordlist.")

        # write ranking of concepts to file
        f = open(self.dataset + '_trebor/paps-'+glm+'.stats','w')
        if 'proto' in self.entries:
            f.write('COGID\tGLID\tCONCEPT\tORIGINS\tPROTO\tREFLEXES\n')
        else:
            f.write('COGID\tGLID\tCONCEPT\tORIGINS\tREFLEXES\n')
        concepts = {}
        for a,b in sorted(paps,key=lambda x:x[1],reverse=True):
            
            a1,a2 = a.split(':')
            a3 = self._id2gl[int(a2)]
            
            try:
                concepts[a3] += [b]
            except:
                concepts[a3] = [b]
            
            # check for number of occurrences
            l = [k for k in self.etd[a] if k != 0]

            # check for proto
            if 'proto' in self.entries:
                proto = self[[k[0] for k in l][0],'proto']
                f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(a1,a2,a3,b,proto,len(l)))
            else:
                f.write('{0}\t{1}\t{2}\t{3}\n'.format(a1,a2,a3,b,len(l)))
        f.close()
        if verbose: print("[i] Wrote stats on paps to file.")

        # write stats on concepts
        f = open(self.dataset+'_trebor/concepts-'+glm+'.stats','w')
        for key in concepts:
            concepts[key] = sum(concepts[key])/len(concepts[key])

        for a,b in sorted(concepts.items(),key=lambda x:x[1],reverse=True):
            f.write('{0}\t{1:.2f}\n'.format(a,b))
        f.close()
        if verbose: print("[i] Wrote stats on concepts to file.")

    def analyze(
            self,
            runs = "default",
            mixed = False,
            verbose = False,
            output_gml = False,
            tar = False,
            usetex = False,
            full_analysis = True,
            plot_dists = False,
            output_plot=False,
            plot_mln = False,
            plot_msn = False,
            **keywords
            ):
        """
        Carry out a full analysis using various parameters.

        Parameters
        ----------
        runs : {str list} (default="default")
            Define a couple of different models to be analyzed. Select between:
            
            * 'default': weighted analysis, using parsimony and weights for
              gains and losses
            * 'topdown': use the traditional approach by
              :evobib:`Nelson-Sathi2011`
            * 'restriction': use the restriction approach
            
            You can also define your own mix of models.

        verbose : bool (default = False)
            If set to c{True}, be verbose when carrying out the analysis.
        usetex : bool (default=True)
            Specify whether you want to use LaTeX to render plots.
        mixed : bool (default=False)
            If set to c{True}, calculate a mixed model by selecting the best
            model for each item separately.
        verbose : bool (default=False)
            Set to c{True} for verbose output.
        output_gml : bool (default=False)
            Set to c{True} in order to output every gain-loss-scenario in
            GML-format.
        full_analysis : bool (default=True)
            Specifies whether a full analysis is carried out or not.
        plot_mln : bool (default=True)
            Select or unselect output plot for the MLN.
        plot_msn : bool (default=False)
            Select or unselect output plot for the MSN.

        """
        
        # set defaults
        defaults = {
                "colorbar" : None, #mpl.cm.jet,
                'threshold':1,
                'fileformat':'pdf',
                'usetex':False,
                'only':[],
                'colormap': None, #mpl.cm.jet
                }

        for key in defaults:
            if key not in keywords:
                keywords[key] = defaults[key]

        # define a default set of runs
        if runs in ['default','weighted']:
            runs = [
                    ('weighted',(3,1)),
                    ('weighted',(5,2)),
                    ('weighted',(2,1)),
                    ('weighted',(3,2)),
                    ('weighted',(1,1)),
                    ]

        elif runs in ['topdown','top-down']:
            runs = [('topdown',2),
                    ('topdown',3),
                    ('topdown',4),
                    ('topdown',5),
                    ]

        elif runs == 'restriction':

            runs = [('restriction',2),
                    ('restriction',3),
                    ('restriction',4),
                    ('restriction',5),
                    ]
        
        # carry out the various analyses
        for mode,params in runs:
            if mode == 'weighted':
                if verbose: print(
                        "[i] Analysing dataset with mode {0} ".format(mode)+\
                                "and ratio {0[0]}:{0[1]}...".format(params)
                                )

                self.get_GLS(
                        mode = mode,
                        ratio = params,
                        verbose = verbose,
                        output_gml = output_gml,
                        tar = tar,
                        output_plot=output_plot
                        )
            elif mode == 'restriction':
                if verbose: print(
                        "[i] Analysing dataset with mode {0} ".format(mode)+\
                                "and restriction {0}...".format(params)
                                )
                
                self.get_GLS(
                        mode = mode,
                        restriction = params,
                        verbose = verbose,
                        output_gml = output_gml,
                        tar = tar,
                        output_plot=output_plot
                        )
            elif mode == 'topdown':
                if verbose: print(
                        "[i] Analysing dataset with mode {0} ".format(mode)+\
                                "and restriction {0}...".format(params)
                                )
                self.get_GLS(
                        mode = mode,
                        restriction = params,
                        verbose = verbose,
                        output_gml = output_gml,
                        tar = tar,
                        output_plot = output_plot
                        )
    
        # calculate the different distributions
        # start by calculating the contemporary distributions
        if verbose: print("[i] Calculating the Contemporary Vocabulary Distributions...")
        self.get_CVSD(verbose=verbose)
        
    
        # now calculate the rest of the distributions
        if verbose: print("[i] Calculating the Ancestral Vocabulary Distributions...")
        modes = list(self.gls.keys())
        for m in modes:
            self.get_AVSD(m,verbose=verbose)

        # calculate mixed model
        if mixed:
            if verbose: print("[i] Calculating the mixed model...")
            self.get_IVSD(
                    verbose=verbose,
                    output_plot=output_plot,
                    output_gml=output_gml,
                    tar=tar
                    )
            modes += ['mixed']

        # compare the distributions using mannwhitneyu
        if verbose: print("[i] Comparing the distributions...")
        
        zp_vsd = []
        for m in modes:
            vsd = sps.mannwhitneyu(
                    self.dists['contemporary'],
                    self.dists[m]
                    )

            zp_vsd.append(vsd)

        # write results to file
        if verbose: print("[i] Writing stats to file.")
        f = open(self.dataset+'_trebor/'+self.dataset+'.stats','w')
        f.write("Mode\tANO\tMNO\tVSD_z\tVSD_p\n")
        for i in range(len(zp_vsd)):
            f.write(
                    '{0}\t{1:.2f}\t{2}\t{3}\n'.format(
                        modes[i],
                        self.stats[modes[i]]['ano'],
                        self.stats[modes[i]]['mno'],
                        '{0[0]}\t{0[1]:.4f}'.format(zp_vsd[i])
                        )
                    )
        f.close()

        # plot the stats if this is defined in the settings
        if plot_dists:

            # specify latex
            mpl.rc('text',usetex=usetex)
                        
            # store distributions in lists
            dists_vsd = [self.dists[m] for m in modes]
            
            # store contemporary dists
            dist_vsd = self.dists['contemporary']
            
            # get the average number of origins
            ano = [self.stats[m]['ano'] for m in modes]

            # create a sorter for the distributions
            sorter = [s[0] for s in sorted(
                zip(range(len(modes)),ano),
                key=lambda x:x[1]   
                )]

            # sort the stuff
            dists_vsd = [dists_vsd[i] for i in sorter]
            modes = [modes[i] for i in sorter]
            mode_strings = [m for m in modes]

            # sort the zp-values
            zp_vsd = [zp_vsd[i] for i in sorter]

            # format the zp-values
            if usetex:

                p_vsd = []
                for i,(z,p) in enumerate(zp_vsd):
                    if p < 0.001:
                        p_vsd.append('p$<${0:.2f}'.format(p))
                    elif p >= 0.05:
                        p_vsd.append(r'\textbf{{p$=${0:.2f}}}'.format(p))
                        
                        # adjust the modes
                        mode_strings[i] = r'\textbf{'+modes[i]+'}'
                    else:
                        p_vsd.append('p$=${0:.2f}'.format(p))
                
            else:
                p_vsd = []
                for z,p in zp_vsd:
                    if p < 0.001:
                        p_vsd.append('p<{0:.2f}'.format(p))
                    elif p >= 0.05:
                        p_vsd.append(r'p={0:.2f}'.format(p))
                    else:
                        p_vsd.append('p={0:.2f}'.format(p))
            
            # create the figure
            fig = plt.figure()

            # create the axis
            ax = fig.add_subplot(111)

            # add the boxplots
            b = ax.boxplot([dist_vsd]+dists_vsd)
            plt.setp(b['medians'],color='black')
            plt.setp(b['whiskers'],color='black')
            plt.setp(b['boxes'],color='black')

            # adjust the yticks
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(18)

            # add the xticks
            plt.xticks(
                    range(1,len(modes)+2),
                    ['']+['{0}\n{1}'.format(
                        m,
                        p
                        ) for m,p in zip(mode_strings,p_vsd)
                        ],
                    size=18
                    )

            # save the figure
            plt.savefig(self.dataset+'_trebor/vsd.pdf')
            plt.clf()
            
            if verbose: print("[i] Plotted the distributions.")
        

        # carry out further analyses if this is specified
        if full_analysis:
            
            # if the mixed model is not chosen
            if not mixed:
                # determine the best model
                p_vsd = [p for z,p in zp_vsd]
                maxP = max(p_vsd)
                glm = modes[p_vsd.index(maxP)]
            else:
                glm = 'mixed'
            
            # set the best model
            self.best_model = glm

            self.get_MLN(
                glm,
                verbose = verbose,
                threshold = keywords['threshold'],
                #colormap = keywords['colormap']
                )

            # check whether plots are chosen
            if plot_mln:
                self.plot_MLN(
                        glm,
                        verbose=verbose,
                        filename=self.dataset+'_trebor/mln-'+glm,
                        threshold = keywords['threshold'],
                        fileformat = keywords['fileformat'],
                        usetex = keywords['usetex'],
                        colormap = keywords['colormap']
                        )
            if plot_msn:
                self.plot_MSN(
                        glm,
                        verbose=verbose,
                        filename=self.dataset+'_trebor/msn-'+glm,
                        fileformat=keywords['fileformat'],
                        threshold = keywords['threshold'],
                        only = keywords['only'],
                        usetex = keywords['usetex'],
                        colormap = keywords['colormap']
                        )

            self.get_PDC(
                    glm,
                    verbose = verbose
                    )

    def plot_MLN(
            self,
            glm = '',
            filename = '',
            fileformat = 'pdf',
            threshold = 1,
            usetex = True,
            colormap = None, #mpl.cm.jet,
            taxon_labels = 'taxon.short_labels',
            verbose = False,
            **keywords
            ):
        """
        Plot the MLN with help of Matplotlib.

        glm : str (default='')
            Identifier for the gain-loss model that is plotted. Defaults to the
            model that had the best scores in terms of probability.
        filename : str (default='')
            If no filename is selected, the filename is identical with the
            dataset.
        fileformat : {'svg','png','jpg','pdf'} (default='pdf')
            Select the format of the output plot.
        threshold : int (default=1)
            Select the threshold for drawing lateral edges.
        verbose : bool (default = False)
            If set to c{True}, be verbose when carrying out the analysis.
        usetex : bool (default=True)
            Specify whether you want to use LaTeX to render plots.
        colormap : {None matplotlib.cm}
            A :py:class:`matplotlib.colormap` instance. If set to c{None}, this
            defaults to :py:class:`~matplotlib.cm.jet`.
        taxon_labels : str (default='taxon.short_labels')
            Specify the taxon labels that should be included in the plot. 

        """
        # check for correct glm
        if not glm and hasattr(self,'best_model'):
            glm = self.best_model
        elif not glm:
            raise ValueError(
            "[i] You should select an appropriate model first."
            )
        
        # check for filename
        if not filename:
            filename = self.dataset

        # if not colormap
        if not colormap:
            colormap = mpl.cm.jet

        # set default, XXX change later
        if 'height' not in keywords:
            keywords['height'] = 7
        if 'width' not in keywords:
            keywords['width'] = 10

        # try to load the configuration file
        try:
            conf = json.load(open(self.dataset+'.json'))
        except:
            conf = {}
        
        # check for 'taxon.labels' in conf
        if taxon_labels in conf: #XXX change later
            tfunc = lambda x:conf[taxon_labels][x]
        else:
            tfunc = lambda x:x

        # get the graph
        graph = self.graph[glm]

        # store in internal and external nodes
        inodes = []
        enodes = []

        # get some data on the taxa
        max_label_len = max([len(tfunc(t)) for t in self.taxa])

        # get colormap for edgeweights
        edge_weights = []
        for nodeA,nodeB,data in graph.edges(data=True):
            if data['label'] == 'horizontal': 
                edge_weights.append(data['weight'])

        
        # determine a colorfunction
        cfunc = np.array(np.linspace(10,256,len(set(edge_weights))),dtype='int')
        lfunc = np.linspace(0.5,8,len(set(edge_weights)))

        # sort the weights
        weights = sorted(set(edge_weights))

        # get the scale for the weights (needed for the line-width)
        scale = 20.0 / max(edge_weights)

        # write colors and scale to graph
        for nA,nB,data in graph.edges(data=True):
            if data['label'] == 'horizontal':
                w = data['weight']
                data['graphics'] = {}
                data['graphics']['fill'] = mpl.colors.rgb2hex(colormap(cfunc[weights.index(w)]))
                data['graphics']['width'] = scale * w

        # get the nodes
        for n,d in graph.nodes(data=True):
            g = d['graphics']
            x = g['x']
            y = g['y']
            h = g['h']
            w = g['w']

            if d['label'] not in self.taxa:
                inodes += [(x,y)]
            else:
                if 'angle' in d['graphics']:
                    r = d['graphics']['angle']
                else:
                    r = 0

                # get the difference between the current label and it's
                # original lenght for formatting output
                ll = max_label_len - len(tfunc(d['label']))

                if usetex:
                    enodes += [(
                        x,
                        y,
                        r'\textbf{'+tfunc(d['label']).replace('_',r'\_')+r'}',
                        r
                        )]
                else:
                    enodes += [(x,y,tfunc(d['label']),r)]
        
        # store vertical and lateral edges
        vedges = []
        ledges = []
        weights = []

        # get the edges
        for a,b,d in graph.edges(data=True):
            
            xA = graph.node[a]['graphics']['x']
            yA = graph.node[a]['graphics']['y']
            xB = graph.node[b]['graphics']['x']
            yB = graph.node[b]['graphics']['y']

            if d['label'] == 'vertical':

                vedges += [(xA,xB,yA,yB)]
            else:
                g = d['graphics']
                f = g['fill']
                w = g['width']
                if d['weight'] < threshold:
                    w = 0.0

                ledges += [(xA,xB,yA,yB,f,w)]

                weights.append(d['weight'])
        
        # usetex
        mpl.rc('text',usetex = usetex)

        # create the figure
        fig = plt.figure(
                facecolor='white',
                figsize = (keywords['width'],keywords['height'])
                    )
        figsp = fig.add_subplot(111)
        
        # create the axis
        ax = plt.axes(frameon=False)
        plt.xticks([0],[''])
        plt.yticks([0],[''])
        
        # set equal axis
        plt.axis('equal')

        # draw the horizontal edges
        for xA,xB,yA,yB,f,w in sorted(ledges,key=lambda x:x[-1]):
            plt.plot(
                    [xA,xB],
                    [yA,yB],
                    '-',
                    color=f,
                    linewidth=float(w) / 3,
                    alpha=0.75
                    )

        # draw the vertical edges
        for xA,xB,yA,yB in vedges:
            plt.plot(
                    [xA,xB],
                    [yA,yB],
                    '-',
                    color='0.0',
                    linewidth=5,
                    )
            plt.plot(
                    [xA,xB],
                    [yA,yB],
                    '-',
                    color='1.0',
                    linewidth=2,
                    )
        # store x,y values for ylim,xlim drawing
        xvals = []
        yvals = []

        # draw the nodes
        for x,y in inodes:
            xvals += [x]
            yvals += [y]

            plt.plot(
                    x,
                    y,
                    'o',
                    markersize=10,
                    color='black',
                    )
            plt.plot(
                    x,
                    y,
                    'o',
                    markersize=6,
                    color='white'
                    )
        
        # draw the leaves
        # store x and y-maxima for ylim, xlim drawing
        for x,y,t,r in enodes:
            
            xvals += [x]
            yvals += [y]

            if r > 270:
                ha = 'right'
            else:
                ha = 'left'
            plt.text(
                    x,
                    y,
                    t,
                    size = '7',
                    verticalalignment='baseline',
                    backgroundcolor='black',
                    horizontalalignment=ha,
                    fontweight = 'bold',
                    color='white',
                    bbox = dict(
                        facecolor='black',
                        boxstyle='square'
                        ),
                    rotation=r,
                    rotation_mode = 'anchor'
                    )

        # add a colorbar
        cax = figsp.imshow(
                [[1,2],[1,2]],
                cmap=colormap,
                visible=False
                )
        cbar = fig.colorbar(
                cax,
                ticks = [
                    1,
                    1.25,
                    1.5,
                    1.75,
                    2
                    ],
                orientation='vertical',
                shrink=0.55
                )
        cbar.set_clim(1.0)
        cbar.set_label('Inferred Links')
        cbar.ax.set_yticklabels(
                [
                    str(min(weights)),
                    '',
                    str(int(max(weights) / 2)),
                    '',
                    str(max(weights))
                    ]
                )
        plt.ylim(min(yvals),max(yvals))
        plt.xlim(min(xvals),max(xvals))
        plt.subplots_adjust(left=0.05,right=0.99,top=0.95,bottom=0.00)
        #fig.axes.get_xaxis().set_visible(False)
        #fig.axes.get_yaxis().set_visible(False)
        plt.axis('off')


        # save the figure
        plt.savefig(filename+'.'+fileformat,bbbox_inches='tight')
        plt.clf()
        if verbose: FileWriteMessage(filename,fileformat).message('written')

        return 

    def plot_MSN(
            self,
            glm = '',
            verbose=False,
            filename='pdf',
            fileformat='pdf',
            threshold = 1,
            only = [],
            usetex = True,
            external_edges = False,
            colormap = None, #mpl.cm.jet
            ):
        """
        Plot the Minimal Spatial Network.

        Parameters
        ----------
        glm : str (default='')
            A string that encodes which model should be plotted.
        filename : str
            The name of the file to which the plot shall be written.
        fileformat : str
            The output format of the plot.
        threshold : int (default=1)
            The threshold for the minimal amount of shared links that shall be
            plotted.
        only : list (default=[])
            The list of taxa whose connections with other taxa should be
            plotted.
        usetex : bool (default=True)
            Specify whether LaTeX shall be used for the plot.

        """
        # check for correct glm
        if not glm and hasattr(self,'best_model'):
            glm = self.best_model
        elif not glm:
            raise ValueError(
            "[i] You should select an appropriate model first."
            )
        
        # check for filename
        if not filename:
            filename = self.dataset

        # check for only
        if not only:
            only = self.taxa

        if not colormap:
            colormap = mpl.cm.jet

        # usetex
        mpl.rc('text',usetex=True)
    
        # redefine taxa and tree for convenience
        taxa,tree = self.taxa,self.tree

        # get the graph
        graph = self.graph[glm]

        # XXX check for coordinates of the taxa, otherwise load them from file and
        # add them to the wordlist XXX add later, we first load it from file
        if 'coords' in self.entries:
            pass
        
        else:
            coords = csv2dict(
                    self.dataset,
                    'coords',
                    dtype=[str,float,float]
                    )

        # check for groups, add functionality for groups in qlc-file later XXX
        if 'group' in self.entries:
            pass
        else:
            groups = dict([(k,v) for k,v in csv2list(self.dataset,'groups')])
        ## check for color, add functionality for colors later XXX
        if 'colors' in self.entries:
            pass
        else:
            colors = dict([(k,v) for k,v in csv2list(self.dataset,'colors')])

        if verbose: LoadDataMessage('coordinates','groups','colors').message('loaded')
        
        # load the rc-file XXX add internal loading later
        try:
            conf = json.load(open(self.dataset+'.json'))
        except:
            pass # XXX add fallback later
        
        if verbose: LoadDataMessage('configuration')
                
        # calculate all resulting edges, using convex hull as
        # approximation 
        geoGraph = nx.Graph()
        
        for nA,nB,d in graph.edges(data=True):
            
            # get the labels
            lA = graph.node[nA]['label']
            lB = graph.node[nB]['label']
            
            # first check, whether edge is horizontal
            if d['label'] == 'horizontal':
                
                # if both labels occur in taxa, it is simple
                if lA in taxa and lB in taxa:
                    try:
                        geoGraph.edge[lA][lB]['weight'] += d['weight']
                    except:
                        geoGraph.add_edge(lA,lB,weight=d['weight'])
                elif not external_edges:
                    # if only one in taxa, we need the convex hull for that node
                    if lA in taxa or lB in taxa:

                        # check which node is in taxa
                        if lA in taxa:
                            this_label = lA
                            other_nodes = tree.getNodeMatchingName(lB).getTipNames()
                            other_label = lB
                        elif lB in taxa:
                            this_label = lB
                            other_nodes = tree.getNodeMatchingName(lA).getTipNames()
                            other_label = lA

                        # get the convex points of others
                        these_coords = [(round(coords[t][0],5),round(coords[t][1],5)) for t in
                                other_nodes]
                        hulls = getConvexHull(these_coords,polygon=False)
    
                        # get the hull with the minimal euclidean distance
                        distances = []
                        for hull in hulls:
                            distances.append(linalg.norm(np.array(hull) - np.array(coords[this_label])))
                        this_hull = hulls[distances.index(min(distances))]
                        other_label = other_nodes[
                                these_coords.index(
                                    (
                                        round(this_hull[0],5),
                                        round(this_hull[1],5)
                                        )
                                    )
                                ]
    
                        # append the edge to the graph
                        try:
                            geoGraph.edge[this_label][other_label]['weight'] += d['weight']
                        except:
                            geoGraph.add_edge(this_label,other_label,weight=d['weight'])
                        
                    #else:
                    #    # get the taxa of a and b
                    #    taxA = tree.getNodeMatchingName(lA).getTipNames()
                    #    taxB = tree.getNodeMatchingName(lB).getTipNames()
    
                    #    # get the convex points
                    #    coordsA = [(round(coords[t][0],5),round(coords[t][1],5)) for t in taxA]
                    #    coordsB = [(round(coords[t][0],5),round(coords[t][1],5)) for t in taxB]
                    #    hullsA = getConvexHull(coordsA,polygon=False)
                    #    hullsB = getConvexHull(coordsB,polygon=False)
    
                    #    # get the closest points
                    #    distances = []
                    #    hulls = []
                    #    for i,hullA in enumerate(hullsA):
                    #        for j,hullB in enumerate(hullsB):
                    #            distances.append(linalg.norm(np.array(hullA)-np.array(hullB)))
                    #            hulls.append((hullA,hullB))
                    #    minHulls = hulls[distances.index(min(distances))]
                    #    
                    #    labelA = taxA[coordsA.index((round(minHulls[0][0],5),round(minHulls[0][1],5)))]
                    #    labelB = taxB[coordsB.index((round(minHulls[1][0],5),round(minHulls[1][1],5)))]
                    #    
                    #    # append the edge to the graph
                    #    try:
                    #        geoGraph.edge[labelA][labelB]['weight'] += d['weight']
                    #    except:
                    #        geoGraph.add_edge(labelA,labelB,weight=d['weight'])

        # get the weights for the lines
        weights = []
        for a,b,d in geoGraph.edges(data=True):
            weights += [d['weight']]
        max_weight = max(weights)
        scale = 256 / max_weight
        sorted_weights = sorted(set(weights))

        # get a color-function
        color_dict = np.array(
                np.linspace(
                    0,
                    256,
                    len(set(weights))
                    ),
                dtype='int'
                )

        # get a line-function
        line_dict = np.linspace(
                0.5,
                conf['linewidth'],
                len(set(weights))
                )

        # scale the weights for line-widths
        linescale = conf['linescale'] / (max_weight-threshold) #XXX
        # XXX apparently not needed?
        
        # determine the maxima of the coordinates
        latitudes = [i[0] for i in coords.values()]
        longitudes = [i[1] for i in coords.values()]

        min_lat,max_lat = min(latitudes),max(latitudes)
        min_lon,max_lon = min(longitudes),max(longitudes)

        # start to initialize the basemap
        fig = plt.figure()
        figsp = fig.add_subplot(111)
        
        # instantiate the basemap
        m = bmp.Basemap(
            llcrnrlon=min_lon + conf['min_lon'],
            llcrnrlat=min_lat + conf['min_lat'],
            urcrnrlon=max_lon + conf['max_lon'],
            urcrnrlat=max_lat + conf['max_lat'],
            resolution=conf['resolution'],
            projection=conf['projection']
            )
        
        # draw first values
        m.drawmapboundary(fill_color=conf['water_color'])
        m.drawcoastlines(color=conf['continent_color'],linewidth=0.5)
        m.drawcountries(color=conf['coastline_color'],linewidth=0.5)
        m.fillcontinents(color=conf['continent_color'],lake_color=conf['water_color'])

        # plot the lines
        for a,b,d in sorted(geoGraph.edges(data=True),key=lambda x:x[2]['weight']):
            
            # don't draw lines beyond threshold
            if d['weight'] < threshold:
                pass
            else:
                if a in coords and b in coords and a in only or b in only:
                    w = d['weight']

                    # retrieve the coords
                    yA,xA = coords[a]
                    yB,xB = coords[b]
                    
                    # get the points on the map
                    xA,yA = m(xA,yA)
                    xB,yB = m(xB,yB)

                    # plot the points
                    plt.plot(
                            [xA,xB],
                            [yA,yB],
                            '-',
                            color=colormap(
                                color_dict[sorted_weights.index(w)]
                                ),
                            alpha = conf['alpha'],
                            linewidth=line_dict[sorted_weights.index(w)],
                            zorder = w + 50
                            )

        # plot the points for the languages
        cell_text = []
        legend_check = []

        # check for taxon.labels in conf
        if 'taxon.labels' in conf:
            tfunc = lambda x:conf['taxon.labels'][x]
        else:
            tfunc = lambda x:x
        if 'groups.labels' in conf:
            gfunc = lambda x:conf['groups.labels'][x]
        else:
            gfunc = lambda x:x

        for i,(taxon,(lng,lat)) in enumerate(sorted(coords.items(),key=lambda x:x[0])):
            
            # retrieve x and y from the map
            x,y = m(lat,lng)
            
            # get the color of the given taxon
            #taxon_color = colors[groups[taxon]]
            
            # get colors from conf
            this_group = groups[taxon]
            taxon_color = conf['groups.colors'][this_group]
            taxon_marker = conf['groups.markers'][this_group]


            # check for legend

            if gfunc(groups[taxon]) in legend_check:
                # plot the marker
                plt.plot(
                    x,
                    y,
                    taxon_marker,
                    markersize = conf['markersize'],
                    color = taxon_color,
                    zorder = max_weight+52,
                    )
            else:
                # plot the marker
                plt.plot(
                    x,
                    y,
                    taxon_marker,
                    markersize = conf['markersize'],
                    color = taxon_color,
                    zorder = max_weight+52,
                    label=gfunc(groups[taxon])
                    )
                legend_check.append(gfunc(groups[taxon]))
            
            # add number to celltext
            if usetex:
                cell_text.append([str(i+1),tfunc(taxon).replace('_',r'\_')])
            else:
                cell_text.append([str(i+1),tfunc(taxon)])

            # plot the text
            # check for darkness of color
            if taxon_color in ['black','gray'] or taxon_color[:3] in ['0.3','0.2','0.1','0.0']:
                text_color = 'white'
            else:
                text_color = 'black'

            plt.text(
                x,
                y,
                str(i+1),
                size = str(int(conf['markersize'] / 2)),
                color = text_color,
                label=taxon,
                horizontalalignment='center',
                fontweight="bold",
                verticalalignment='center',
                zorder=max_weight+55
                )

        # add a colorbar
        cax = figsp.imshow(
                [[1,2],[1,2]],
                visible=False,
                cmap=colormap
                )
        cbar = fig.colorbar(
                cax,
                ticks = [
                    1,
                    1.25,
                    1.5,
                    1.75,
                    2
                    ],
                orientation='vertical',
                shrink=0.55
                )
        cbar.set_clim(1.0)
        cbar.set_label('Inferred Links')
        cbar.ax.set_yticklabels(
                [
                    str(min(weights)),
                    '',
                    str(int(max(weights) / 2)),
                    '',
                    str(max(weights))
                    ]
                )

        # add the legend
        this_table = plt.table(
                cellText = cell_text,
                colWidths = conf['table.column.width'],
                loc = conf['table.location'],
                )
        this_table.auto_set_font_size(False)
        this_table.set_fontsize(conf['table.text.size'])

        # adjust the table
        for line in this_table._cells:
            this_table._cells[line]._text._horizontalalignment = 'left'
            this_table._cells[line]._text._fontproperties.set_weight('bold')
            this_table._cells[line]._text.set_color(conf['table.text.color'])
            this_table._cells[line].set_height(conf['table.cell.height'])
            #this_table._cells[line]._text._fontproperties.set_size(conf['table.text.size'])
            this_table._cells[line].set_linewidth(0.0)
            this_table._cells[line].set_color(conf['table.cell.color'])
        
        this_table.set_zorder(100)
        
        plt.legend(
                loc=conf['legend.location'],
                numpoints=1,
                prop={
                    'size':conf['legend.size'],
                    'weight':'bold'
                    }
                )

        plt.subplots_adjust(left=0.02,right=0.98,top=1.0,bottom=0.00)

        plt.savefig(filename+'.'+fileformat)
        plt.clf()
        if verbose: FileWriteMessage(filename,fileformat).message('written')
        return

    def plot_concepts(
            self,
            concept,
            cogA,
            cogB,
            labels = {1:'1',2:'2',3:'3',4:'4'},
            tcolor = {
                1:'white',
                2:'black',
                3:'0.5',
                4:'0.1'
                },
            verbose=False,
            filename='pdf',
            fileformat='pdf',
            threshold = 1,
            usetex = True
            ):
        """
        Plot the Minimal Spatial Network.

        Parameters
        ----------
        glm : str
            A string that encodes which model should be plotted.
        filename : str
            The name of the file to which the plot shall be written.
        fileformat : str
            The output format of the plot.
        threshold : int (default=1)
            The threshold for the minimal amount of shared links that shall be
            plotted.
        only : list (default=[])
            The list of taxa whose connections with other taxa should be
            plotted.
        usetex : bool (default=True)
            Specify whether LaTeX shall be used for the plot.

        """
        # usetex
        mpl.rc('text',usetex=True)
    
        # redefine taxa and tree for convenience
        taxa,tree = self.taxa,self.tree

        # XXX check for coordinates of the taxa, otherwise load them from file and
        # add them to the wordlist XXX add later, we first load it from file
        if 'coords' in self.entries:
            pass
        
        else:
            coords = csv2dict(
                    self.dataset,
                    'coords',
                    dtype=[str,float,float]
                    )

        # check for groups, add functionality for groups in qlc-file later XXX
        if 'group' in self.entries:
            pass
        else:
            groups = dict([(k,v) for k,v in csv2list(self.dataset,'groups')])
        # check for color, add functionality for colors later XXX
        #if 'colors' in self.entries:
        #    pass
        #else:
        #    colors = dict([(k,v) for k,v in csv2list(self.dataset,'colors')])

        if verbose: LoadDataMessage('coordinates','groups','colors').message('loaded')
        
        # load the rc-file XXX add internal loading later
        try:
            conf = json.load(open(self.dataset+'.json'))
        except:
            pass # XXX add fallback later
        
        if verbose: LoadDataMessage('configuration')
                
        # get the paps
        these_taxa = {}
        for taxon in taxa:

            # get the dictionary and the entry
            try:
                cogs = self.get_dict(col=taxon,entry='pap')[concept]
            except:
                cogs = []
            
            # check for identical cogs and assign them to the 4 categories
            if cogA in cogs and cogB in cogs:
                these_taxa[taxon] = 3
            elif cogA in cogs and cogB not in cogs:
                these_taxa[taxon] = 1
            elif cogA not in cogs and cogB in cogs:
                these_taxa[taxon] = 2
            else:
                these_taxa[taxon] = 4

        # determine the maxima of the coordinates
        latitudes = [i[0] for i in coords.values()]
        longitudes = [i[1] for i in coords.values()]

        min_lat,max_lat = min(latitudes),max(latitudes)
        min_lon,max_lon = min(longitudes),max(longitudes)

        # start to initialize the basemap
        fig = plt.figure()
        figsp = fig.add_subplot(111)
        
        # instantiate the basemap
        m = bmp.Basemap(
            llcrnrlon=min_lon + conf['min_lon'],
            llcrnrlat=min_lat + conf['min_lat'],
            urcrnrlon=max_lon + conf['max_lon'],
            urcrnrlat=max_lat + conf['max_lat'],
            resolution=conf['resolution'],
            projection=conf['projection']
            )
        
        # draw first values
        m.drawmapboundary(fill_color=conf['water_color'])
        m.drawcoastlines(color=conf['continent_color'],linewidth=0.5)
        m.drawcountries(color=conf['coastline_color'],linewidth=0.5)
        m.fillcontinents(color=conf['continent_color'],lake_color=conf['water_color'])

        # plot the points for the languages
        cell_text = []
        legend_check = []
        for i,(taxon,(lng,lat)) in enumerate(sorted(coords.items(),key=lambda x:x[0])):
            
            # retrieve x and y from the map
            x,y = m(lat,lng)
            
            # get the color of the given taxon
            #taxon_color = colors[groups[taxon]]
            
            if these_taxa[taxon] == 4:
                marker = '*'
            else:
                marker = 's'

            # check for legend
            if labels[these_taxa[taxon]] in legend_check:
                # check for forth marker
                # plot the marker
                plt.plot(
                    x,
                    y,
                    marker,
                    markersize = conf['markersize'],
                    color = tcolor[these_taxa[taxon]],
                    #zorder = 50,
                    )
            else:
                # plot the marker
                plt.plot(
                    x,
                    y,
                    marker,
                    markersize = conf['markersize'],
                    color = tcolor[these_taxa[taxon]],
                    #zorder = 52,
                    label=labels[these_taxa[taxon]]
                    )
                legend_check.append(labels[these_taxa[taxon]])
            
            # add number to celltext
            if usetex:
                cell_text.append([str(i+1),taxon.replace('_',r'\_')])
            else:
                cell_text.append([str(i+1),taxon])

            # plot the text
            if tcolor[these_taxa[taxon]] == 'black':
                textcolor = 'white'
            else:
                textcolor='black'
            
            plt.text(
                x,
                y,
                str(i+1),
                size = str(int(conf['markersize'] / 2)),
                label=taxon,
                color = textcolor,
                horizontalalignment='center',
                verticalalignment='center',
                )

        this_table = plt.table(
                cellText = cell_text,
                colWidths = conf['table.column.width'],
                loc = conf['table.location'],
                )

        # adjust the table
        for line in this_table._cells:
            this_table._cells[line]._text._horizontalalignment = 'left'
            this_table._cells[line]._text._fontproperties.set_weight('bold')
            this_table._cells[line]._text.set_color(conf['table.text.color'])
            this_table._cells[line].set_height(conf['table.cell.height'])
            this_table._cells[line]._text._fontproperties.set_size(conf['table.text.size'])
            this_table._cells[line].set_linewidth(0.0)
            this_table._cells[line].set_color(conf['table.cell.color'])
        
        this_table.set_zorder(100)
        
        plt.legend(
                loc=conf['legend.location'],
                numpoints=1,
                prop={
                    'size':conf['legend.size'],
                    'weight':'bold'
                    }
                )

        plt.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05)

        plt.savefig(filename+'.'+fileformat)
        plt.clf()
        if verbose: FileWriteMessage(filename,fileformat).message('written')
        return
            


#!depr    def get_weighted_GLS_multi(
#!depr            self,
#!depr            pap,
#!depr            ratio = (1,1),
#!depr            verbose = False
#!depr            ):
#!depr        """
#!depr        Calculate a weighted gain-loss-scenario (WGLS) for a given PAP.
#!depr
#!depr        In contrast to the other modes, bifurcating trees are allowed.
#!depr        """
#!depr        
#!depr        # make a dictionary that stores the scenario
#!depr        d = {}
#!depr        
#!depr        # get the list of nodes that are not missing
#!depr        taxa,paps = [],[]
#!depr        for i,taxon in enumerate(self.taxa):
#!depr            if pap[i] != -1:
#!depr                taxa += [taxon]
#!depr                paps += [pap[i]]
#!depr
#!depr        # get the subtree of all taxa
#!depr        tree = self.tree.getSubTree(taxa)
#!depr
#!depr        # get the subtree containing all taxa that have positive paps
#!depr        tree = tree.lowestCommonAncestor(
#!depr                [
#!depr                    self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1
#!depr                    ]
#!depr                )
#!depr
#!depr        if verbose: print("[i] Subtree is {0}.".format(str(tree)))
#!depr
#!depr        # assign the basic (starting) values to the dictionary
#!depr        nodes = [t.Name for t in tree.tips()]
#!depr
#!depr        if verbose: print("[i] Nodes are {0}.".format(','.join(nodes)))
#!depr
#!depr        # get the first state of all nodes and store the state in the
#!depr        # dictionary. note that we start from two distinct scenarios: one
#!depr        # assuming single origin at the root where all present states in the
#!depr        # leave are treated as retentions, and one assuming multiple origins,
#!depr        # where all present states in the leaves are treated as origins
#!depr        for node in nodes:
#!depr            idx = taxa.index(node)
#!depr            if paps[idx] >= 1:
#!depr                state = 1
#!depr            else:
#!depr                state = 0
#!depr            d[node] = [(state,[])]
#!depr
#!depr        # return simple scenario, if the group is single-origin
#!depr        if sum([d[node][0][0] for node in nodes]) == len(nodes):
#!depr            return [(tree.Name,1)]
#!depr
#!depr        # order the internal nodes according to the number of their leaves
#!depr        ordered_nodes = sorted(
#!depr                tree.nontips()+[tree],key=lambda x:len(x.tips())
#!depr                )
#!depr
#!depr        # calculate the general restriction value (maximal weight). This is roughly
#!depr        # spoken simply the minimal value of either all events being counted as
#!depr        # origins (case 1) or assuming origin of the character at the root and
#!depr        # counting all leaves that lost the character as single loss events (case
#!depr        # 2). In case two, the first gain of the character has to be added
#!depr        # additionally
#!depr        maxWeight = min(paps.count(1) * ratio[0], paps.count(0) * ratio[1] + ratio[0])
#!depr
#!depr        # join the nodes successively
#!depr        for i,node in enumerate(ordered_nodes):
#!depr            if verbose: print(node.Name)
#!depr            
#!depr            # when dealing with multifurcating trees, we have to store all
#!depr            # possible scenarios, i.e. we need to store the crossproduct of all
#!depr            # scenarios
#!depr
#!depr            # get the names of the children of the nodes
#!depr            names = [x.Name for x in node.Children]
#!depr
#!depr            # get the nodes with their states from the dictionary
#!depr            tmp_nodes = [d[x.Name] for x in node.Children]
#!depr
#!depr            # get the cross-product of the stuff
#!depr            crossp = itertools.product(*tmp_nodes)
#!depr            
#!depr            newNodes = []
#!depr            
#!depr            # combine the histories of the items if all have the same value,
#!depr            # therefore, we first get the states in a simple list
#!depr            for cross in crossp:
#!depr                states = [x[0] for x in cross]
#!depr                stories = [x[1] for x in cross]
#!depr                states_sum = sum(states)
#!depr                states_len = len(states)
#!depr                
#!depr                # combine the histories
#!depr                new_stories = []
#!depr                for x in stories:
#!depr                    new_stories += x
#!depr
#!depr                if states_sum == states_len or states_sum == 0:
#!depr
#!depr                    # add the histories to the queue only if their weight is
#!depr                    # less or equal to the maxWeight
#!depr                    gl = [k[1] for k in new_stories]+[x for x in [states[0]] if x == 1]
#!depr
#!depr                    weight = gl.count(1) * ratio[0] + gl.count(0) * ratio[1]
#!depr                    if weight <= maxWeight:
#!depr                        newNodes.append((states[0],new_stories))
#!depr
#!depr                # if the states are not identical, we check for both scenarios
#!depr                else:
#!depr                    # first scenario (tmpA) assumes origin, that is, for each node
#!depr                    # that has a 1, we add an origin to new_stories, same is
#!depr                    # for loss scenario (tmpB)
#!depr                    tmpA = [x for x in new_stories]
#!depr                    tmpB = [x for x in new_stories]
#!depr                    for c,state in enumerate(states):
#!depr                        if state == 1:
#!depr                            tmpA += [(names[c],1)]
#!depr                        else:
#!depr                            tmpB += [(names[c],0)]
#!depr
#!depr                    # get the vectors to make it easier to retrieve the number
#!depr                    # of losses and gains
#!depr                    glA = [k[1] for k in tmpA]
#!depr                    glB = [k[1] for k in tmpB] + [1] # don't forget adding 1 origin
#!depr
#!depr                    # check the gain-loss scores
#!depr                    weightA = glA.count(1) * ratio[0] + glA.count(0) * ratio[1]
#!depr                    weightB = glB.count(1) * ratio[0] + glB.count(0) * ratio[1]
#!depr
#!depr                    newNodeA = (0,tmpA)
#!depr                    newNodeB = (1,tmpB)
#!depr
#!depr                    if 1 in [k[1] for k in tmpB]:
#!depr                        noB = True
#!depr                    else:
#!depr                        noB = False
#!depr
#!depr                    if weightA <= maxWeight:
#!depr                        newNodes += [newNodeA]
#!depr                    if weightB <= maxWeight and not noB:
#!depr                        newNodes += [newNodeB]
#!depr                        
#!depr                        # check whether an additional gain is inferred on either of
#!depr                        # the two possible paths. 
#!depr                        # XXX reduce later, this is not efficient XXX
#!depr                        #if nodeB[0] == 1 and 1 in [k[1] for k in tmpA]:
#!depr                        #    noA = True
#!depr                        #else:
#!depr                        #    noA = False
#!depr                        #if nodeA[0] == 1 and 1 in [k[1] for k in tmpB]:
#!depr                        #    noB = True
#!depr                        #else:
#!depr                        #    noB = False
#!depr              
#!depr                        ## if the weights are above the theoretical maximum, discard
#!depr                        ## the solution
#!depr                        #if weightA <= maxWeight and not noA:
#!depr                        #    newNodes += [newNodeA]
#!depr                        #if weightB <= maxWeight and not noB:
#!depr                        #    newNodes += [newNodeB]
#!depr
#!depr                d[node.Name] = newNodes
#!depr                if verbose: print("nodelen",len(d[node.Name]))
#!depr        
#!depr        # try to find the best scenario by counting the ratio of gains and losses.
#!depr        # the key idea here is to reduce the number of possible scenarios according
#!depr        # to a given criterion. We choose the criterion of minimal changes as a
#!depr        # first criterion to reduce the possibilities, i.e. we weight both gains
#!depr        # and losses by 1 and select only those scenarios where gains and losses
#!depr        # sum up to a minimal number of gains and losses. This pre-selection of
#!depr        # scenarios can be further reduced by weighting gains and losses
#!depr        # differently. So in a second stage we choose only those scenarios where
#!depr        # there is a minimal amount of gains. 
#!depr        
#!depr        if verbose: print(len(d[tree.Name]))
#!depr
#!depr        # convert the specific format of the d[tree.Name] to simple format
#!depr        gls_list = []
#!depr        for first,last in d[tree.Name]:
#!depr            if first == 1:
#!depr                gls_list.append([(tree.Name,first)]+last)
#!depr            else:
#!depr                gls_list.append(last)
#!depr
#!depr        # the tracer stores all scores
#!depr        tracer = []
#!depr
#!depr        for i,line in enumerate(gls_list):
#!depr            
#!depr            # calculate gains and losses
#!depr            gains = sum([1 for x in line if x[1] == 1])
#!depr            losses = sum([1 for x in line if x[1] == 0])
#!depr
#!depr            # calculate the score
#!depr            score = ratio[0] * gains + ratio[1] * losses
#!depr
#!depr            # append it to the tracer
#!depr            tracer.append(score)
#!depr        
#!depr        # get the minimum score
#!depr        minScore = min(tracer)
#!depr
#!depr        # return the minimal indices, sort them according to the number of
#!depr        # gains inferred, thereby pushing gains to the root, similar to
#!depr        # Mirkin's (2003) suggestion
#!depr        return sorted(
#!depr                [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore],
#!depr                key = lambda x:sum([i[1] for i in x])
#!depr                )[0]
#!depr
#!depr
#!depr
#!depr    def get_restricted_GLS(
#!depr            self,
#!depr            pap,
#!depr            restriction = 4,
#!depr            verbose = False
#!depr            ):
#!depr        """
#!depr        Calculate a restricted gain-loss-scenario (RGLS) for a given PAP.
#!depr    
#!depr        """
#!depr
#!depr        # make a dictionary that stores the scenario
#!depr        d = {}
#!depr
#!depr        # get the subtree containing all taxa that have positive paps
#!depr        tree = self.tree.lowestCommonAncestor(
#!depr                [self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1]
#!depr                )
#!depr
#!depr        # assign the basic (starting) values to the dictionary
#!depr        nodes = [x.Name for x in tree.tips()]
#!depr
#!depr        # get the first state of all nodes and store the state in the dictionary.
#!depr        # note that we start from two distinct scenarios: one assuming single
#!depr        # origin at the root, where all present states in the leave are treated as
#!depr        # retentions, and one assuming multiple origins, where all present states
#!depr        # in the leaves are treated as origins
#!depr        for node in nodes:
#!depr            idx = self.taxa.index(node)
#!depr            if pap[idx] >= 1:
#!depr                state = 1
#!depr            else:
#!depr                state = 0
#!depr            d[node] = [(state,[])]
#!depr
#!depr        # return simple scenario if the group is single-origin
#!depr        if sum([d[node][0][0] for node in nodes]) == len(nodes):
#!depr            return [(tree.Name,1)]
#!depr
#!depr        # order the internal nodes according to the number of their leaves
#!depr        ordered_nodes = sorted(
#!depr                tree.nontips()+[tree],key=lambda x:len(x.tips())
#!depr                )
#!depr
#!depr        # join the nodes successively
#!depr        for i,node in enumerate(ordered_nodes):
#!depr            if verbose: print(node)
#!depr            
#!depr            # get the name of the children of the nodes
#!depr            nameA,nameB = [x.Name for x in node.Children]
#!depr
#!depr            # get the nodes with their states from the dictionary
#!depr            nodesA,nodesB = [d[x.Name] for x in node.Children]
#!depr
#!depr            # there might be alternative states, therefore it is important to
#!depr            # iterate over all possible paths
#!depr            newNodes = []
#!depr            for nodeA in nodesA:
#!depr                for nodeB in nodesB:
#!depr                    # if the nodes have the same state, the state is assigned to
#!depr                    # the most recent common ancestor node
#!depr                    if nodeA[0] == nodeB[0]:
#!depr
#!depr                        # combine the rest of the histories of the items
#!depr                        tmp = nodeA[1] + nodeB[1]
#!depr
#!depr                        # append only if tmp is in concordance with maximal_gains
#!depr                        if len([k for k in tmp if k[1] == 1]) + nodeA[0] <= restriction:
#!depr                            newNodes.append((nodeA[0],tmp))
#!depr
#!depr                    # if the nodes have different states, we go on with two
#!depr                    # distinct scenarios
#!depr                    else:
#!depr                        
#!depr                        # first scenario assumes retention of nodeA
#!depr                        tmpA = nodeA[1] + nodeB[1] + [(nameA,nodeA[0])]
#!depr
#!depr                        # second scenario assumes retention of nodeB
#!depr                        tmpB = nodeA[1] + nodeB[1] + [(nameB,nodeB[0])]
#!depr
#!depr                        newNodeA = nodeB[0],tmpA
#!depr                        newNodeB = nodeA[0],tmpB
#!depr                        
#!depr                        # check whether one of the solutions is already above the
#!depr                        # maximum / best scenario with respect to the gain-loss
#!depr                        # weights as a criterion
#!depr
#!depr                        # first, calculate gains and losses
#!depr                        gainsA = nodeB[0]+sum([k[1] for k in tmpA])
#!depr                        gainsB = nodeA[0]+sum([k[1] for k in tmpB])
#!depr
#!depr                        # check whether an additional gain is inferred on either of
#!depr                        # the two possible paths. 
#!depr                        # XXX reduce later, this is not efficient XXX
#!depr                        if nodeB[0] == 1 and 1 in [k[1] for k in tmpA]:
#!depr                            noA = True
#!depr                        else:
#!depr                            noA = False
#!depr                        if nodeA[0] == 1 and 1 in [k[1] for k in tmpB]:
#!depr                            noB = True
#!depr                        else:
#!depr                            noB = False
#!depr                                            
#!depr                        # if the gains are about the theoretical maximum, discard
#!depr                        # the solution
#!depr                        if gainsA <= restriction and not noA:
#!depr                            newNodes += [newNodeA]
#!depr                        if gainsB <= restriction and not noB:
#!depr                            newNodes += [newNodeB]
#!depr
#!depr                d[node.Name] = newNodes
#!depr        
#!depr        # try to find the best scenario by counting the ratio of gains and losses.
#!depr        # the key idea here is to reduce the number of possible scenarios according
#!depr        # to a given criterion. We choose the criterion of minimal changes as a
#!depr        # first criterion to reduce the possibilities, i.e. we weight both gains
#!depr        # and losses by 1 and select only those scenarios where gains and losses
#!depr        # sum up to a minimal number of gains and losses. This pre-selection of
#!depr        # scenarios can be further reduced by weighting gains and losses
#!depr        # differently. So in a second stage we choose only those scenarios where
#!depr        # there is a minimal amount of gains. 
#!depr        
#!depr        if verbose: print(len(d[tree.Name]))
#!depr
#!depr        # convert the specific format of the d[tree.Name] to simple format
#!depr        gls_list = []
#!depr        for first,last in d[tree.Name]:
#!depr            if first == 1:
#!depr                gls_list.append([(tree.Name,first)]+last)
#!depr            else:
#!depr                gls_list.append(last)
#!depr
#!depr        # the tracer stores all scores
#!depr        tracer = []
#!depr
#!depr        for i,line in enumerate(gls_list):
#!depr            
#!depr            # calculate gains and losses
#!depr            gains = sum([1 for x in line if x[1] == 1])
#!depr            losses = sum([1 for x in line if x[1] == 0])
#!depr
#!depr            # calculate the score
#!depr            score = gains + losses
#!depr
#!depr            # append it to the tracer
#!depr            tracer.append(score)
#!depr        
#!depr        # get the minimum score
#!depr        minScore = min(tracer)
#!depr
#!depr        # push gains down to the root as suggested by Mirkin 2003
#!depr        minimal_gains = [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore]
#!depr        
#!depr        best_scenario = None
#!depr        old_length_of_tips = len(self.taxa) + 1
#!depr
#!depr        for i,line in enumerate(minimal_gains):
#!depr            
#!depr            # calculate number of tips for the gains of a given scenario
#!depr            new_length_of_tips = 0
#!depr            for taxon,state in line:
#!depr                if state == 1:
#!depr                    new_length_of_tips += len(
#!depr                            self.tree.getNodeMatchingName(taxon).getTipNames()
#!depr                            )
#!depr            if new_length_of_tips < old_length_of_tips:
#!depr                old_length_of_tips = new_length_of_tips
#!depr                best_scenario = i
#!depr
#!depr        return minimal_gains[best_scenario]

#!depr    def get_weighted_GLS(
#!depr            self,
#!depr            pap,
#!depr            ratio = (1,1),
#!depr            verbose = False
#!depr            ):
#!depr        """
#!depr        Calculate a weighted gain-loss-scenario (WGLS) for a given PAP.
#!depr        """
#!depr        
#!depr        # make a dictionary that stores the scenario
#!depr        d = {}
#!depr        
#!depr        # get the list of nodes that are not missing
#!depr        taxa,paps = [],[]
#!depr        for i,taxon in enumerate(self.taxa):
#!depr            if pap[i] != -1:
#!depr                taxa += [taxon]
#!depr                paps += [pap[i]]
#!depr
#!depr        # get the subtree of all taxa
#!depr        tree = self.tree.getSubTree(taxa)
#!depr
#!depr        # get the subtree containing all taxa that have positive paps
#!depr        tree = tree.lowestCommonAncestor(
#!depr                [
#!depr                    self.taxa[i] for i in range(len(self.taxa)) if pap[i] >= 1
#!depr                    ]
#!depr                )
#!depr
#!depr        if verbose: print("[i] Subtree is {0}.".format(str(tree)))
#!depr
#!depr        # assign the basic (starting) values to the dictionary
#!depr        nodes = [t.Name for t in tree.tips()]
#!depr
#!depr        if verbose: print("[i] Nodes are {0}.".format(','.join(nodes)))
#!depr
#!depr        # get the first state of all nodes and store the state in the
#!depr        # dictionary. note that we start from two distinct scenarios: one
#!depr        # assuming single origin at the root where all present states in the
#!depr        # leave are treated as retentions, and one assuming multiple origins,
#!depr        # where all present states in the leaves are treated as origins
#!depr        for node in nodes:
#!depr            idx = taxa.index(node)
#!depr            if paps[idx] >= 1:
#!depr                state = 1
#!depr            else:
#!depr                state = 0
#!depr            d[node] = [(state,[])]
#!depr
#!depr        # return simple scenario, if the group is single-origin
#!depr        if sum([d[node][0][0] for node in nodes]) == len(nodes):
#!depr            return [(tree.Name,1)]
#!depr
#!depr        # order the internal nodes according to the number of their leaves
#!depr        ordered_nodes = sorted(
#!depr                tree.nontips()+[tree],key=lambda x:len(x.tips())
#!depr                )
#!depr
#!depr        # calculate the general restriction value (maximal weight). This is roughly
#!depr        # spoken simply the minimal value of either all events being counted as
#!depr        # origins (case 1) or assuming origin of the character at the root and
#!depr        # counting all leaves that lost the character as single loss events (case
#!depr        # 2). In case two, the first gain of the character has to be added
#!depr        # additionally
#!depr        maxWeight = min(paps.count(1) * ratio[0], paps.count(0) * ratio[1] + ratio[0])
#!depr
#!depr        # join the nodes successively
#!depr        for i,node in enumerate(ordered_nodes):
#!depr            if verbose: print(node.Name)
#!depr            
#!depr            # get the name of the children of the nodes
#!depr            nameA,nameB = [x.Name for x in node.Children]
#!depr
#!depr            # get the nodes with their states from the dictionary
#!depr            nodesA,nodesB = [d[x.Name] for x in node.Children]
#!depr
#!depr            # there might be alternative states, therefore it is important to
#!depr            # iterate over all possible paths
#!depr            newNodes = []
#!depr            for nodeA in nodesA:
#!depr                for nodeB in nodesB:
#!depr                    # if the nodes have the same state, the state is assigned to
#!depr                    # the most recent common ancestor node
#!depr                    if nodeA[0] == nodeB[0]:
#!depr
#!depr                        # combine the rest of the histories of the items
#!depr                        tmp = nodeA[1] + nodeB[1]
#!depr
#!depr                        # add them to the queue only if their weight is less or
#!depr                        # equal to the maxWeight
#!depr                        gl = [k[1] for k in tmp]+[x for x in [nodeA[0]] if x == 1]
#!depr                        weight = gl.count(1) * ratio[0] + gl.count(0) * ratio[1]
#!depr
#!depr                        if weight <= maxWeight:
#!depr                            newNodes.append((nodeA[0],tmp))
#!depr
#!depr                    # if the nodes have different states, we go on with two
#!depr                    # distinct scenarios
#!depr                    else:
#!depr                        
#!depr                        # first scenario assumes retention of nodeB
#!depr                        tmpA = nodeA[1] + nodeB[1] + [(nameA,nodeA[0])]
#!depr
#!depr                        # second scenario assumes retention of nodeA
#!depr                        tmpB = nodeA[1] + nodeB[1] + [(nameB,nodeB[0])]
#!depr
#!depr                        # get the vectors in order to make it easier to retrieve
#!depr                        # the number of losses and gains
#!depr                        glA = [k[1] for k in tmpA] + [x for x in [nodeB[0]] if x == 1]
#!depr                        glB = [k[1] for k in tmpB] + [x for x in [nodeA[0]] if x == 1]
#!depr
#!depr                        # check the gain-loss scores 
#!depr                        weightA = glA.count(1) * ratio[0] + glA.count(0) * ratio[1]
#!depr                        weightB = glB.count(1) * ratio[0] + glB.count(0) * ratio[1]
#!depr                        
#!depr                        # check whether an additional gain is inferred on either of
#!depr                        # the two possible paths. 
#!depr                        # XXX reduce later, this is not efficient XXX
#!depr                        if nodeB[0] == 1 and 1 in [k[1] for k in tmpA]:
#!depr                            noA = True
#!depr                        else:
#!depr                            noA = False
#!depr                        if nodeA[0] == 1 and 1 in [k[1] for k in tmpB]:
#!depr                            noB = True
#!depr                        else:
#!depr                            noB = False
#!depr
#!depr                        newNodeA = nodeB[0],tmpA
#!depr                        newNodeB = nodeA[0],tmpB
#!depr                                            
#!depr                        # if the weights are above the theoretical maximum, discard
#!depr                        # the solution
#!depr                        if weightA <= maxWeight and not noA:
#!depr                            newNodes += [newNodeA]
#!depr                        if weightB <= maxWeight and not noB:
#!depr                            newNodes += [newNodeB]
#!depr
#!depr                d[node.Name] = newNodes
#!depr                if verbose: print("nodelen",len(d[node.Name]))
#!depr        
#!depr        # try to find the best scenario by counting the ratio of gains and losses.
#!depr        # the key idea here is to reduce the number of possible scenarios according
#!depr        # to a given criterion. We choose the criterion of minimal changes as a
#!depr        # first criterion to reduce the possibilities, i.e. we weight both gains
#!depr        # and losses by 1 and select only those scenarios where gains and losses
#!depr        # sum up to a minimal number of gains and losses. This pre-selection of
#!depr        # scenarios can be further reduced by weighting gains and losses
#!depr        # differently. So in a second stage we choose only those scenarios where
#!depr        # there is a minimal amount of gains. 
#!depr        
#!depr        if verbose: print(len(d[tree.Name]))
#!depr
#!depr        # convert the specific format of the d[tree.Name] to simple format
#!depr        gls_list = []
#!depr        for first,last in d[tree.Name]:
#!depr            if first == 1:
#!depr                gls_list.append([(tree.Name,first)]+last)
#!depr            else:
#!depr                gls_list.append(last)
#!depr
#!depr        # the tracer stores all scores
#!depr        tracer = []
#!depr
#!depr        for i,line in enumerate(gls_list):
#!depr            
#!depr            # calculate gains and losses
#!depr            gains = sum([1 for x in line if x[1] == 1])
#!depr            losses = sum([1 for x in line if x[1] == 0])
#!depr
#!depr            # calculate the score
#!depr            score = ratio[0] * gains + ratio[1] * losses
#!depr
#!depr            # append it to the tracer
#!depr            tracer.append(score)
#!depr        
#!depr        # get the minimum score
#!depr        minScore = min(tracer)
#!depr
#!depr        # return the minimal indices, sort them according to the number of
#!depr        # gains inferred, thereby pushing gains to the root, similar to
#!depr        # Mirkin's (2003) suggestion
#!depr        return sorted(
#!depr                [gls_list[i] for i in range(len(tracer)) if tracer[i] == minScore],
#!depr                key = lambda x:sum([i[1] for i in x])
#!depr                )[0]
    

