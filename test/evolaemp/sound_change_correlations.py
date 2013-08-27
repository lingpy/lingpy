from numpy import *
from lingpy2 import *

# switch namespace to evolaemp
from lingpy2.data.names.evolaemp import *
from lingpy2.align.sca import get_consensus
from lingpy2.thirdparty import cogent as cg
from ete2 import Tree, TextFace

import math

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

#load long language names
f = open('data/asjp/world_longnames.txt','r')
rl = f.readlines()
f.close()
longnames = array([x.strip() for x in rl])
nameToID = dict({(longnames[i],i) for i in range(0,len(longnames))})

#load long language names
f = open('data/asjp/world_names.txt','r')
rl = f.readlines()
f.close()
names = array([x.strip() for x in rl])

guideTree = cg.LoadTree("data/asjp/world-NWPV.nwk")
#convert guideTree node names to integers as expected by Lingpy MSA
for leaf in guideTree.tips():
    leaf.Name = str(nameToID[leaf.Name])   

langs = range(1603,1670) #iranian languages
langs = range(1426,1459) #germanic languages
langs = range(1461,1602) #indic languages
langs = range(1690,1707) #slavic languages
langs = range(1671,1690) #romance languages
phylName = "IE"
familyName = "ROMANCE"
familyGuideTree = subGuideTree(guideTree,langs)
familyNameTable = [longnames[lang] for lang in langs]
nameTable = [names[lang] for lang in langs]
#printTree(familyGuideTree,0,names=familyNameTable)

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
    lexstat = LexStat(lexdict,model=internal_asjp,merge_vowels=False)
    lexstat.get_scorer()
    lexstat.cluster(method='sca',threshold=0.55,verbose=True)
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
        cognateLangs = [int(lexdict[IDList[0]][0]) - langs[0] for IDList in etym_dict[cognateID] if IDList != 0]
        if len(cognateLangs) > 1:  #cognate sets of size 1 are useless
            cognateGuideTree = subGuideTree(familyGuideTree,cognateLangs)       
            #print("\nAligning cognate " + str(cognateID) + ":")
            #print "  cognate langs = " + str(cognateLangs)
            #printTree(cognateGuideTree,0,names=[germanicNameTable[lang] for lang in cognateLangs])
            cognateNameTable = [nameTable[lang] for lang in cognateLangs]
            multi = MSA("./cognate" + str(cognateID) + ".msq",merge_vowels=True,unique_seqs=False)
            tree_mtx = convert.newick.nwk2guidetree(str(cognateGuideTree))
            multi.prog_align(model=sca,gop=-4,scale=0.9,guide_tree=tree_mtx)
            #print(multi)
            #old version of call had taxa=[str(i) for i in range(len(langs))]
            cons = get_consensus(multi, cognateGuideTree, gaps=True, classes=False)
            print("Reconstructed proto-" + familyName + " word for concept " + str(conceptID - 3) + ":\t" + cons)  
             
            #PRINT OUT RECONSTRUCTION STEPS IN A TREE VISUALIZATION
            outputTree = Tree()
            outputTree.add_face(TextFace(str("".join(cognateGuideTree.reconstructed))), column=0, position = "branch-right")
            def copyChildrenIntoOutput(treeNode, outputTreeNode):
                for child in treeNode.Children:
                    outputChild = outputTreeNode.add_child()
                    if child.isTip():
                        outputChild.name = str("".join(child.reconstructed)) + " (" + cognateNameTable[int(child.Name)] + ")"
                    else:
                        outputChild.add_face(TextFace(str("".join(child.reconstructed))), column=0, position = "branch-right")
                        copyChildrenIntoOutput(child, outputChild) 
            copyChildrenIntoOutput(cognateGuideTree,outputTree)
            outputTree.render("output/" + phylName + "/" + familyName + "/" + str(conceptID - 3) + "." + cons.replace("-","") +".png")
            
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
printTree(familyGuideTree,0,names=familyNameTable, field="recon_changes")