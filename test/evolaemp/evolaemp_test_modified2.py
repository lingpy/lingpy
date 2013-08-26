from numpy import *
from lingpy2 import *

# switch namespace to evolaemp
from lingpy2.data.names.evolaemp import *
from lingpy2.align.sca import get_consensus
from lingpy2.thirdparty import cogent as cg

import math

#READING IN THE ASJP DATA

f = open('data/asjp/sounds41.txt','r')
sounds = array([x.strip() for x  in f.readlines()])
f.close()

# reading in 'asjpMatrix.txt' into a numpy array
mfile = open("data/asjp/asjpMatrix.txt","r")
asjpRaw = mfile.readlines()
mfile.close()

asjpMatrix = array([x.strip().split('\t') for x in asjpRaw])


# restricting asjpMatrix to the languages in world_names.txt
f = open('data/asjp/world_names.txt','r')
rl = f.readlines()
f.close()

names = array([x.strip() for x in rl])

asjpMatrix = array([x for x in asjpMatrix if x[0] in names])

internal_asjp = Model("asjp_el")

#TEST 1: LINGPY-BASED LANGUAGE DISTANCE MEASURE
print("\nTest 1: language distance measure based on pairwise alignment")
print("----------------------------------------------------------")

# Lingpy-based distance between row l1 and row l2 in mtx
def ldistLingpy(l1,l2,mtx=asjpMatrix):
    l1List = mtx[l1,4:]
    l2List = mtx[l2,4:]
    pairs = []
    wordOffsets = []
    for wordID in range(0,len(l1List)):
      if l1List[wordID] != '0' and l2List[wordID] != '0':
        wordOffsets.append(len(pairs))
        for entry1 in l1List[wordID].split('-'):
          for entry2 in l2List[wordID].split('-'):
            pairs.append((entry1,entry2))
    wordOffsets.append(len(pairs))
    #align all word pairs in parallel
    pair = Pairwise(pairs,merge_vowels=False)
    pair.align(distance=True,model=asjp)
    #collect the lowest distance values for each wordID (i.e. concept)
    distValues = [min(alignment[2] for alignment in pair.alignments[wordOffsets[i]:wordOffsets[i+1]]) for i in range(0,len(wordOffsets)-1)]
    #simply compute the average distance value
    averageDist = sum(distValues) / len(distValues)
    print("  ldistLingpy(" + asjpMatrix[l1][0] + "," + asjpMatrix[l2][0] + ") = " + str(averageDist))

ldistLingpy(1,2)
ldistLingpy(3,10)
ldistLingpy(5,6)
ldistLingpy(120,150)
ldistLingpy(235,1810)

#TEST 2: EXTRACTING PHONEME REPLACEMENT COUNTS USING MSA
print("\nTest 2: Extracting Phoneme Replacement Counts Using MSA")
print("----------------------------------------------------------")

lexdict = {}
lexdict[0] = ["ID", "concept", "ipa", "doculect"]
ID = 1

langs = range(1426,1459) #germanic languages

replacementOccurrenceTable = dict((asjpMatrix[taxon1,0],{}) for taxon1 in langs)
for taxon1 in langs:
  for taxon2 in langs:
    replacementOccurrenceTable[asjpMatrix[taxon1,0]][asjpMatrix[taxon2,0]] = {} 

#create a dictionary for cognate detection
for langID in langs:
  entries = asjpMatrix[langID,39] #entry for "mountain"
  for entry in entries.split('-'):
    if not entry == '0':
      lexdict[ID] = [langID, "mountain", entry, asjpMatrix[langID,0]]
      ID += 1

#cluster words into cognate sets
lexstat = LexStat(lexdict,model=internal_asjp,merge_vowels=False)
lexstat.get_scorer()
lexstat.cluster(method='sca',threshold=0.5,verbose=True)
etym_dict = lexstat.get_etymdict(ref='scaid', entry='', loans=False)

for cognateID in etym_dict.keys():
    entry_msq_file = open("cognate" + str(cognateID) + ".msq", 'w')
    entry_msq_file.write("ASJP database\n")
    entry_msq_file.write("Cognate " + str(cognateID) + " for Germanic languages\n")
    for IDList in etym_dict[cognateID]:
      if (IDList != 0):
        [langID, word, entry, langName] = lexdict[IDList[0]][:4]
        entry_msq_file.write(langName + "\t" + entry + "\n")
    entry_msq_file.close()
    print("Aligning cognate " + str(cognateID) + ":\n")
    multi = MSA("./cognate" + str(cognateID) + ".msq",merge_vowels=False)
    multi.prog_align(model=sca,gop=-2,scale=0.7)
    print(multi)
    #collect the sound replacements in this cognate set
    cognateSize = len(multi.seqs)
    if cognateSize > 1:
      mtx = multi.alm_matrix
      length = len(mtx[0])
      for i in range(0,cognateSize):
        for j in range(0,cognateSize):
          for k in range(0,length):
            #print(multi.taxa[i] + "\t->\t" + multi.taxa[j] + ":\t" + mtx[i][k] + "/" + mtx[j][k])
            occurrenceDict = replacementOccurrenceTable[multi.taxa[i]][multi.taxa[j]]
            #print(occurrenceDict)
            if mtx[i][k] not in occurrenceDict:
              occurrenceDict[mtx[i][k]] = {}
            if mtx[j][k] not in occurrenceDict[mtx[i][k]]:
              occurrenceDict[mtx[i][k]][mtx[j][k]] = 0.0
            occurrenceDict[mtx[i][k]][mtx[j][k]] += 1.0
            
#normalize the numbers in the occurrence dictionary
for taxon1 in replacementOccurrenceTable.keys():
  for taxon2 in replacementOccurrenceTable[taxon1].keys():
    for phon1 in replacementOccurrenceTable[taxon1][taxon2].keys():
      entries = replacementOccurrenceTable[taxon1][taxon2][phon1]
      entrySum = sum(list(entries.values()))
      for phon2 in entries.keys():
        print("  " + taxon1 + "->" + taxon2 + ", " + phon1 + "->" + phon2 + ": " + str(entries[phon2]) + "/" + str(entrySum))
        entries[phon2] = entries[phon2] / entrySum
        
#TEST 3: EXTRACTION OF SUB-GUIDETREES WITH POINTERS TO THE ORIGINAL
print("\nTest 3: Extraction of Sub-Guidetrees with Pointers to the Original")
print("-------------------------------------------------------------------------")
def selectNodes(tree, selIndices):
    selNodes = []
    for leaf in tree.tips():
        if int(leaf.Name) in selIndices:
            selNodes.append(leaf)
    return selNodes

def treePath(node):
    path = [node]
    while not node.isRoot():
        node = node.Parent
        path.insert(0,node)
    return path

def constructSubtree(paths,index,curNode,indexMap):
    #create a map [node -> all paths containing that node at index position]
    partition = {}
    for node in {path[index] for path in paths}:
        partition[node] = [path for path in paths if path[index] == node]
    #partition = {(node,[path for path in paths if path[index] == node]) for node in {path[index] for path in paths}}
    if len(partition) == 1:
        #no split, we simply go on to the next index in paths
        constructSubtree(paths,index + 1,curNode,indexMap)
    else:
        #split according to the partition, creating a new node where necessary
        for node in partition.keys():
            if len(partition[node]) == 1:
                #we have arrived at a leaf (or a unary branch above it), copy the leaf
                newLeafName = str(indexMap[int(partition[node][0][-1].Name)]) 
                newLeaf = cg.tree.TreeNode(Name=newLeafName)
                newLeaf.orig = partition[node][0][-1]
                curNode.Children.append(newLeaf)
                newLeaf.Parent = curNode
            else:               
                newNode = cg.tree.TreeNode()
                newNode.orig = node
                curNode.Children.append(newNode)
                newNode.Parent = curNode
                constructSubtree(partition[node],index + 1,newNode,indexMap)         

def subGuideTree(tree,selIndices):
    selNodes = selectNodes(tree,selIndices)
    indexMap = dict(zip(selIndices,range(0,len(selIndices))))
    paths = [treePath(node) for node in selNodes]
    #print str(paths)
    subtree = cg.tree.TreeNode()
    subtree.orig = tree.root()
    constructSubtree(paths,1,subtree,indexMap)
    return subtree

def printTree(node, depth, names=None, field=None, func=None):
    name = node.Name
    if name == None:
        name = "no name"
    else:
        if names != None:
            name = name + ": " + names[int(name)]
    if field != None and hasattr(node,field):
        if func != None:
            name = name + " " + func(getattr(node,field))
        else:
            name = name + " " + str(getattr(node,field))
    name = depth * "  " + name;
    print name;
    for child in node.Children:
        printTree(child, depth + 1, names, field, func)

def printTreeWithOrigPointers(node, depth, names=None):
    name = node.Name
    if name == None:
        name = "no name"
    else:
        if names != None:
            name = names[int(name)]
    name = depth * "  " + name;
    orig = node
    while hasattr(orig, "orig"):
        name += " -> " + str(orig.orig)
        orig = orig.orig
    print name;
    for child in node.Children:
        printTreeWithOrigPointers(child, depth + 1, names)

print "Original Tree: " + "(((0,1),(2,(3,4))),(5,6));"
print "\nSelected subnodes: [0,2,3,5,6]"
tree = cg.LoadTree(treestring="(((0,1),(2,(3,4))),(5,6));")
subtree = subGuideTree(tree,[0,2,3,5,6])
print "Subtree: " + str(subtree)
print "With pointers:"
printTreeWithOrigPointers(subtree,0)
print "\nSubselected nodes: [1,2,3]"
subsubtree = subGuideTree(subtree,[1,2,3])
print "Subsubtree: " + str(subsubtree)
print "With pointers:"
printTreeWithOrigPointers(subsubtree,0)

#TEST 4: USER-DEFINED GUIDE TREES FOR MSA
print("\nTest 4: User-Defined Guide Trees for MSA")
print("----------------------------------------------------------")

multi1 = MSA("data/asjp/test-alignment.msq",merge_vowels=False)
print("\nAlignment with clustering-based guide tree:")
multi1.prog_align(model=sca,gop=-2,scale=0.7)
print("Inferred guide tree: " + str(multi1.tree_matrix))
print(multi1)
print("\nAlignment with custom guide tree:")
multi2 = MSA("data/asjp/test-alignment.msq",merge_vowels=False)
tree_mtx = convert.newick.nwk2guidetree("(((0,2),(1,3)),(5,4));")
multi2.prog_align(model=sca,gop=-2,scale=0.7,guide_tree=tree_mtx)
print("User-defined guide tree: " + str(multi2.tree_matrix))
print(multi2)
print("\nWrite partial alignments and sizes into the tree:")
tree = cg.LoadTree(treestring="(((0,2),(1,3)),(5,4));")
for node in tree.postorder():
    if node.isTip():
        node.alignment = [multi2.alm_matrix[int(node.Name)]]
        node.size = 1
    else:
        node.alignment = node.Children[0].alignment + node.Children[1].alignment
        node.size = node.Children[0].size + node.Children[1].size
for node in tree.postorder():
    print node.Name + ": (" + str(node.size) + ") " + str(node.alignment)
print("\nCompute phoneme distribution at each position of the alignment:")
for node in tree.postorder():
    node.distribution = []
for i in range (0,len(multi2.alm_matrix[0])):
    for node in tree.postorder():
        if node.isTip():
            node.distribution.append({multi2.alm_matrix[int(node.Name)][i] : 1.0})
        else:
            node.distribution.append({})
            child1 = node.Children[0]
            child2 = node.Children[1]
            for phoneme in set(child1.distribution[i].keys()) | set(child2.distribution[i].keys()):
                value = 0.0
                if phoneme in child1.distribution[i].keys():
                    value += child1.size * child1.distribution[i][phoneme]
                if phoneme in child2.distribution[i].keys():
                    value += child2.size * child2.distribution[i][phoneme]
                value /= node.size
                node.distribution[i][phoneme] = value
for node in tree.postorder():
    print node.Name + ": " + str(node.distribution)
print("\nReconstruct word forms at inner nodes by simplistic criteria:")
for node in tree.postorder():
    node.reconstructed = []
for i in range (0,len(multi2.alm_matrix[0])):
    for node in tree.postorder():
        dist = node.distribution[i]
        maxValue = max(dist.values())
        maxKeys = [key for key in dist.keys() if dist[key]==maxValue]
        if len(maxKeys) == 1 or node.isRoot():
            node.reconstructed.append(maxKeys[0])
        else:         
            parentDist = node.Parent.distribution[i]
            maxKey = max(maxKeys, key=(lambda key: parentDist[key]))
            node.reconstructed.append(maxKey)
for node in tree.postorder():
    print node.Name + ": " + "".join(node.reconstructed)
print("\nDetermining and counting sound changes at the edges of the guide tree:")
for node in tree.postorder():
    node.recon_changes = {}
    if not node.isRoot():
        for i in range (0,len(node.reconstructed)):
            if node.reconstructed[i] != node.Parent.reconstructed[i]:
                change = (node.Parent.reconstructed[i], node.reconstructed[i])
                if change not in node.recon_changes.keys():
                    node.recon_changes[change] = 1
                else:
                    node.recon_changes[change] += 1
for node in tree.postorder():
    print node.Name + ": " + str(node.recon_changes)

#TEST 5: COMBINING SUB-GUIDETREE SELECTION AND RECONSTRUCTION
print("\nTest 5: Combining Sub-Guidetree Selection and Reconstruction")
print("--------------------------------------------------------------------")

#load long language names
f = open('data/asjp/world_longnames.txt','r')
rl = f.readlines()
f.close()
longnames = array([x.strip() for x in rl])
nameToID = dict({(longnames[i],i) for i in range(0,len(longnames))})

guideTree = cg.LoadTree("data/asjp/world-NWPV.nwk")
#convert guideTree node names to integers as expected by Lingpy MSA
for leaf in guideTree.tips():
    leaf.Name = str(nameToID[leaf.Name])   

langs = range(1603,1670) #iranian languages
langs = range(1671,1690) #romance languages
langs = range(1690,1707) #slavic languages
langs = range(1461,1602) #indic languages
langs = range(1426,1459) #germanic languages
germanicGuideTree = subGuideTree(guideTree,langs)
germanicNameTable = [longnames[lang] for lang in langs]
printTree(germanicGuideTree,0,names=germanicNameTable)

for conceptID in range(4,44):
    #repeat the cognate detection process to generate cognate sets
    lexdict = {}
    lexdict[0] = ["ID", "concept", "ipa", "doculect"]
    ID = 1
    #create a dictionary for cognate detection
    for langID in langs:
        entries = asjpMatrix[langID,conceptID] #originally: 39 for "mountain"
        for entry in entries.split('-'):
            if not entry == '0':
                lexdict[ID] = [langID, "concept" + str(conceptID), entry, asjpMatrix[langID,0]]
                ID += 1
    
    #cluster words into cognate sets
    lexstat = LexStat(lexdict,model=internal_asjp,merge_vowels=True)
    lexstat.get_scorer()
    lexstat.cluster(method='sca',threshold=0.55,verbose=True)
    etym_dict = lexstat.get_etymdict(ref='scaid', entry='', loans=False)
    
    print("etym_dict.keys() = " + str(etym_dict.keys()))
    for cognateID in etym_dict.keys():
        entry_msq_file = open("cognate" + str(cognateID) + ".msq", 'w')
        entry_msq_file.write("ASJP database\n")
        entry_msq_file.write("Cognate " + str(cognateID) + " for Germanic languages\n")
        for IDList in etym_dict[cognateID]:
            if (IDList != 0):
                [langID, word, entry, langName] = lexdict[IDList[0]][:4]
                entry_msq_file.write(langName + "\t" + entry + "\n")
        entry_msq_file.close()
        cognateLangs = [int(lexdict[IDList[0]][0]) - langs[0] for IDList in etym_dict[cognateID] if IDList != 0]
        if len(cognateLangs) > 1:  #cognate sets of size 1 are useless
            cognateGuideTree = subGuideTree(germanicGuideTree,cognateLangs)       
            print("\nAligning cognate " + str(cognateID) + ":")
            print "  cognate langs = " + str(cognateLangs)
            printTree(cognateGuideTree,0,names=[germanicNameTable[lang] for lang in cognateLangs])
            multi = MSA("./cognate" + str(cognateID) + ".msq",merge_vowels=True,unique_seqs=False)
            tree_mtx = convert.newick.nwk2guidetree(str(cognateGuideTree))
            multi.prog_align(model=sca,gop=-4,scale=0.9,guide_tree=tree_mtx)
            print(multi)
            cons = get_consensus(multi, cognateGuideTree, gaps=True, taxa=[str(i) for i in range(len(langs))], classes=False)
            print("Reconstructed proto word for concept " + str(conceptID - 3) + ":\t" + cons)   
            #print("\nDetermining and counting sound changes at the edges of the guide tree, and cascading them to the supertrees:")
            for node in cognateGuideTree.postorder():
                if not hasattr(node, "recon_changes"):
                    node.recon_changes = {}
                if not node.isRoot():
                    for i in range (0,len(node.reconstructed)):
                        #if node.reconstructed[i] != node.Parent.reconstructed[i]:
                            change = (node.Parent.reconstructed[i], node.reconstructed[i])
                            origNode = node
                            while origNode != None:
                                if not hasattr(origNode, "recon_changes"):
                                    origNode.recon_changes = {}
                                if change not in origNode.recon_changes.keys():
                                    origNode.recon_changes[change] = 1
                                else:
                                    origNode.recon_changes[change] += 1
                                if hasattr(origNode, "orig"):
                                    origNode = origNode.orig
                                else:
                                    origNode = None
            #for node in cognateGuideTree.postorder():
            #    print str(node.Name)
            #    print str(node.Name) + ": " + str(node.recon_changes)
            #printTree(cognateGuideTree,0,names=[germanicNameTable[lang] for lang in cognateLangs], field="recon_changes")

print("\nSupertree with collected sound changes at the edges:")
printTree(germanicGuideTree,0,names=germanicNameTable, field="recon_changes")