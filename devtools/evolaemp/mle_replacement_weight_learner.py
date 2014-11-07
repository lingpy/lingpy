from numpy import *
from lingpy import *
import os

# switch namespace to evolaemp
from lingpy.data.names.evolaemp import *
from lingpy.align.sca import get_consensus
from lingpy.thirdparty import cogent as cg
from ete2 import Tree, TextFace

import math

graphicalOutput = False
cognateDetection = False

method = "strengthenClearWinner"

def ensure_dir(f):
    try:
        os.makedirs(f)
    except OSError:
        pass
    #d = os.path.dirname(f)
    #if not os.path.exists(d):
    #    os.makedirs(d)

def printTree(node, depth, names=None, field=None, func=None):
    name = node.Name
    if name == None:
        name = "no name"
    else:
        if names != None and name != "root":
            name = name + ": " + names[int(name)]
    if field != None and hasattr(node,field):
        if func != None:
            name = name + " " + func(getattr(node,field))
        else:
            name = name + " " + str(getattr(node,field))
    else:
        name = name + "NO FIELD " + field + "!"
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

#f = open('data/asjp/sounds41.txt','r')
#sounds = array([x.strip() for x  in f.readlines()]).tolist()
#f.close()

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
longNameToID = dict({(longnames[i],i) for i in range(0,len(longnames))})

#load long language names
f = open('data/asjp/world_names.txt','r')
rl = f.readlines()
f.close()
names = array([x.strip() for x in rl])
nameToID = dict({(names[i],i) for i in range(0,len(names))})



guideTree = cg.LoadTree("data/asjp/world-NWPV.nwk")
#convert guideTree node names to integers as expected by Lingpy MSA
for leaf in guideTree.tips():
    leaf.Name = str(longNameToID[leaf.Name]) 
    
iteration = 1
numIterations = 100

while iteration <= numIterations:
    if iteration == 1:
        #mfile = open("replacement-weights.txt","r")
        #sounds = array(mfile.readline().strip().split("\t"))
        #repWeightsRaw = mfile.readlines()
        #mfile.close()        
        #repWeights = array([x.strip().split('\t') for x in repWeightsRaw])
        
        sounds = array(['!','3','4','5','7','8','C','E','G','L','N','S','T','X','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','-'])
        repWeights = array([[-log(0.2 / (len(sounds) - 1)) for phon1 in sounds] for phon2 in sounds])
        for i in range(len(sounds)):
            repWeights[i,i] = -log(0.8) 
    else:    
        #load replacement weights for reconstruction
        mfile = open("replacement-weights" + str(iteration -1) + ".txt","r")
        sounds = array(mfile.readline().strip().split("\t"))
        repWeightsRaw = mfile.readlines()
        mfile.close()
        
        repWeights = array([x.strip().split('\t') for x in repWeightsRaw])
        
    if iteration % 2 == 0:
        method = "strengthenClearWinner"
    else:
        method = "strengthenUnclearWinner"
        
    symbolToSoundID = dict([(sounds[i],i) for i in range(len(sounds))])
    def rep_weights(phon1, phon2):
        if len(phon1) == 1 and len(phon2) == 1:
            return double(repWeights[symbolToSoundID[phon1],symbolToSoundID[phon2]])
        return 10
    
    global_replacement_probabilities = dict()
    for phon1 in sounds:
        global_replacement_probabilities[phon1] = dict()
        for phon2 in sounds:
            global_replacement_probabilities[phon1][phon2] = exp(-rep_weights(phon1,phon2))
    
    globalParsimony = 0  
    globalLikelihood = 0
    
    probCorrectionsPlus = dict()
    for sound in sounds:
        probCorrectionsPlus[sound] = dict((sound,0.0) for sound in sounds)
    probCorrectionsMinus = dict()
    for sound in sounds:
        probCorrectionsMinus[sound] = dict((sound,0.0) for sound in sounds)
    
    for familyName in unique(asjpMatrix[:,2]):
        langs = where(asjpMatrix[:,2] == familyName)[0].tolist()
        phylName = str(asjpMatrix[langs[0],1])
        #print phylName + "." + str(familyName)
        if len(langs) == 1: 
            continue
            #print "only one language, skipped"
    
        familyParsimony = 0.0
        familyLikelihood = 0.0
        
        ensure_dir("output/" + phylName + "/" + familyName)
        ensure_dir("cognates/" + phylName + "/" + familyName)
    
    #langs = range(1603,1670) #iranian languages
    #langs = range(1426,1459) #germanic languages
    #langs = range(1461,1602) #indic languages
    #langs = range(1690,1707) #slavic languages
    #langs = range(1671,1690) #romance languages
    #langs = range(1707,1727) #finnic languages
    #langs = range(1734,1743) #mongolic languages
    #langs = range(1743,1765) #tungusic languages
    #langs = range(1765,1821) #turkic languages
    #langs = range(1291,1347) #semitic languages
    
    #phylName = "IE"
        #print langs
        familyGuideTree = subGuideTree(guideTree,langs)
        familyNameTable = [longnames[lang] for lang in langs]
        nameTable = [names[lang] for lang in langs]
        #printTree(familyGuideTree,0,names=familyNameTable)
        
        for conceptID in range(4,44):
            if cognateDetection:
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
                
                newCognateID = 1    
                for cognateID in etym_dict.keys():
                    entry_msq_file = open("cognates/" + phylName + "/" + familyName + "/" + str(conceptID - 3) + "." + str(newCognateID) + ".msq", 'w')
                    entry_msq_file.write("ASJP database\n")
                    entry_msq_file.write("Cognate " + str(cognateID) + " for Germanic languages\n")
                    for IDList in etym_dict[cognateID]:
                        if (IDList != 0):
                            [langID, word, entry, langName] = lexdict[IDList[0]][:4]
                            entry_msq_file.write(langName + "\t" + entry + "\n")
                    entry_msq_file.close()
                    newCognateID += 1
            
            numCognates = 0
            path = "cognates/" + phylName + "/" + familyName      
            for fname in os.listdir(path):
                if fname.startswith(str(conceptID - 3) + "."):
                    numCognates += 1
            
            for newCognateID in range(1,numCognates + 1):
                multi = MSA("cognates/" + phylName + "/" + familyName + "/" + str(conceptID - 3) + "." + str(newCognateID) + ".msq",merge_vowels=False,unique_seqs=False)
                cognateLangs = [nameToID[taxon] - langs[0] for taxon in multi.taxa]
                #cognateLangs = [int(lexdict[IDList[0]][0]) - langs[0] for IDList in etym_dict[cognateID] if IDList != 0]
                #print "cognate set " + str(cognateID) + " - cognate langs: " + str(cognateLangs)
                if len(cognateLangs) > 1:  #cognate sets of size 1 are useless
                    cognateGuideTree = subGuideTree(familyGuideTree,cognateLangs)       
                    #print("\nAligning cognate " + str(cognateID) + ":")
                    #print "  cognate langs = " + str(cognateLangs)
                    #printTree(cognateGuideTree,0,names=[germanicNameTable[lang] for lang in cognateLangs])
                    cognateNameTable = [nameTable[lang] for lang in cognateLangs]
                    tree_mtx = convert.newick.nwk2guidetree(str(cognateGuideTree))
                    multi.prog_align(model=internal_asjp,gop=-4,scale=0.9,guide_tree=tree_mtx)
                    #print(multi)
                    
                    cons = get_consensus(multi, cognateGuideTree, recon_alg="sankoff_parsimony", gaps=True, classes=False, rep_weights = rep_weights, local = "gap")
                    
                    #collect correction estimates based on the assumption that the reconstruction is correct
                    for node in cognateGuideTree.postorder():
                        if not node.isTip():
                            for i in range(len(node.reconstructed)):
                                #OLD VERSION, this lead to locally suboptimal values because it built on the most parsimonious reconstruction!
                                #minValue = node.sankoffTable[i][node.reconstructed[i]]
                                #optimalChar = node.reconstructed[i]
                                minValue = min(node.sankoffTable[i].values())
                                optimalChar = [key for key in node.sankoffTable[i].keys() if node.sankoffTable[i][key]==minValue][0]
                                if len(optimalChar) == 1:
                                    if len(node.sankoffTable[i].keys()) == 1:
                                        secondPlaceValue = 0
                                    else:
                                        secondPlaceValue = 65000 #simulating Integer.maxInt   
                                        for phon in node.sankoffTable[i].keys():
                                            if phon != optimalChar:
                                                if node.sankoffTable[i][phon] < secondPlaceValue:
                                                    secondPlaceValue = node.sankoffTable[i][phon]
                                    #THIS METHOD ATTEMPTS TO ENLARGE THE DISTANCE OF THE WINNER TO THE REST IN CASE OF UNCLEAR OUTCOMES
                                    if method == "strengthenUnclearWinner":
                                        scalingFactor = 0.0
                                        if secondPlaceValue != 0:
                                            scalingFactor = minValue / secondPlaceValue
                                        if len(node.Children[0].reconstructed[i]) == 1: probCorrectionsPlus[optimalChar][node.Children[0].reconstructed[i]] += 1.0 * scalingFactor
                                        if len(node.Children[1].reconstructed[i]) == 1: probCorrectionsPlus[optimalChar][node.Children[1].reconstructed[i]] += 1.0 * scalingFactor
                                        parsimonySumWrong = sum(node.sankoffTable[i].values()) - node.sankoffTable[i][optimalChar]
                                        for phon in node.sankoffTable[i].keys():
                                            if len(phon) == 1 and phon != optimalChar:
                                                if len(node.Children[0].reconstructed[i]) == 1: probCorrectionsMinus[phon][node.Children[0].reconstructed[i]] += (node.sankoffTable[i][phon] / parsimonySumWrong) * scalingFactor
                                                if len(node.Children[1].reconstructed[i]) == 1: probCorrectionsMinus[phon][node.Children[1].reconstructed[i]] += (node.sankoffTable[i][phon] / parsimonySumWrong) * scalingFactor
                                    #THIS METHOD FURTHER STRENGTHENS THE WINNER (AND WEAKENS THE LOSERS) IN CASE OF A CLEAR OUTCOME
                                    elif method == "strengthenClearWinner":
                                        scalingFactor = 1.0
                                        if secondPlaceValue != 0:
                                            scalingFactor = 1.0 - (minValue / secondPlaceValue)
                                        if len(node.Children[0].reconstructed[i]) == 1: probCorrectionsPlus[optimalChar][node.Children[0].reconstructed[i]] += 1.0 * scalingFactor
                                        if len(node.Children[1].reconstructed[i]) == 1: probCorrectionsPlus[optimalChar][node.Children[1].reconstructed[i]] += 1.0 * scalingFactor
                                        parsimonyDifferenceSum = sum([node.sankoffTable[i][phon] - minValue for phon in node.sankoffTable[i].keys()])
                                        #parsimonyDifferenceSum = 0.0 #comment out this line to activate collection of negative evidence
                                        if parsimonyDifferenceSum != 0.0:
                                            for phon in node.sankoffTable[i].keys():
                                                if len(phon) == 1 and phon != optimalChar:
                                                    if len(node.Children[0].reconstructed[i]) == 1: probCorrectionsMinus[phon][node.Children[0].reconstructed[i]] += ((node.sankoffTable[i][phon] - minValue) / parsimonyDifferenceSum) * scalingFactor
                                                    if len(node.Children[1].reconstructed[i]) == 1: probCorrectionsMinus[phon][node.Children[1].reconstructed[i]] += ((node.sankoffTable[i][phon] - minValue) / parsimonyDifferenceSum) * scalingFactor
    
                    #aggregate the parsimony value
                    cognateParsimony = 0.0
                    for i in range(len(cognateGuideTree.reconstructed)):
                        cognateParsimony += min(cognateGuideTree.sankoffTable[i].values())
                    familyParsimony += cognateParsimony
                    
                    #aggregate the likelihood value
                    cognateLikelihood = 0.0
                    for node in cognateGuideTree.postorder():
                        if not node.isTip():
                            for i in range(len(node.reconstructed)):
                                char = node.reconstructed[i]
                                char1 = node.Children[0].reconstructed[i]
                                char2 = node.Children[1].reconstructed[i]
                                if len(char) == 1 and len(char1) == 1 and len(char2) == 1:
                                    cognateLikelihood += rep_weights(char,char1) + rep_weights(char,char2)
                    familyLikelihood += cognateLikelihood
                        
                    #print("Reconstructed proto-" + familyName + " word for concept " + str(conceptID - 3) + ":\t" + cons + "\twith average parsimony " + str(cognateParsimony / len(cognateLangs)))
                     
                    #PRINT OUT RECONSTRUCTION STEPS IN A TREE VISUALIZATION
                    if graphicalOutput:
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
                                (seg1, seg2) = (str(node.Parent.reconstructed[i]), str(node.reconstructed[i]))
                                origNode = node
                                while origNode != None:
                                    if not hasattr(origNode, "recon_changes"):
                                        origNode.recon_changes = {}
                                    if seg1 not in origNode.recon_changes:
                                        origNode.recon_changes[seg1] = {}
                                    if seg2 not in origNode.recon_changes[seg1]:
                                        origNode.recon_changes[seg1][seg2] = 0.0
                                    origNode.recon_changes[seg1][seg2] += 1.0
                                    #if change not in origNode.recon_changes.keys():
                                    #    origNode.recon_changes[change] = 1
                                    #else:
                                    #    origNode.recon_changes[change] += 1
                                    if hasattr(origNode, "orig"):
                                        origNode = origNode.orig
                                    else:
                                        origNode = None
                    #for node in cognateGuideTree.postorder():
                    #    print str(node.Name)
                    #    print str(node.Name) + ": " + str(node.recon_changes)
                    #printTree(cognateGuideTree,0,names=[germanicNameTable[lang] for lang in cognateLangs], field="recon_changes")
        #print("Total parsimony while deriving proto-" + familyName + ": " + str(familyParsimony))  
        globalParsimony += familyParsimony  
        globalLikelihood += familyLikelihood
    print("Iteration " + str(iteration) + " complete - parsimony " + str(globalParsimony) + ", likelihood " + str(globalLikelihood)) 
    
    global_replacements = dict()
    for phon1 in sounds:
        global_replacements[phon1] = dict()
        for phon2 in sounds:
            global_replacements[phon1][phon2] = 0
    
    for node in guideTree.postorder():
      if hasattr(node, "recon_changes"):
          for phon1 in node.recon_changes.keys():
            entries = node.recon_changes[phon1]
            entrySum = sum(entries.values())
            for phon2 in entries.keys():
              if len(phon1) == 1 and len(phon2) == 1: #avoid problems with "87" and the like
                global_replacements[phon1][phon2] += entries[phon2]
                #print("  " + taxon1 + "->" + taxon2 + ", " + phon1 + "->" + phon2 + ": " + str(entries[phon2]) + "/" + str(entrySum))
                entries[phon2] = entries[phon2] / entrySum
              
    global_replacement_probabilities = dict()
    for phon1 in sounds:
        global_replacement_probabilities[phon1] = dict()
        entrySum = sum(global_replacements[phon1].values())
        entrySum += len([phon2 for phon2 in sounds if global_replacements[phon1][phon2] == 0]) * 0.0001
        for phon2 in sounds:
            if entrySum == 0 or global_replacements[phon1][phon2] == 0:
                global_replacement_probabilities[phon1][phon2] = 0.0001
            else:
                global_replacement_probabilities[phon1][phon2] = global_replacements[phon1][phon2] / entrySum
    
    for phon1 in sounds:
        for phon2 in sounds:
            if global_replacement_probabilities[phon1][phon2] != 0.0001:
                print phon1 + phon2 + ": " + str(- log(global_replacement_probabilities[phon1][phon2]))
        print
    
#     exponent = 1.0/(5 * iteration)
#     #re-estimation of replacement probabilities 
#     for phon1 in sounds:
#         for phon2 in sounds:
#             if probCorrectionsPlus[phon1][phon2] != 0 or probCorrectionsMinus[phon1][phon2] != 0:
#                 posEvidence = probCorrectionsPlus[phon1][phon2]
#                 if (posEvidence < 0.1): posEvidence = 0.1
#                 negEvidence = probCorrectionsMinus[phon1][phon2]
#                 if (negEvidence < 0.1): negEvidence = 0.1
#                 correctionFactor = (posEvidence / negEvidence) ** exponent
#                 #correctionFactor = 1.0
#                 #if probCorrectionsMinus[phon1][phon2] > 1:
#                 #    correctionFactor *= 1/(log(probCorrectionsMinus[phon1][phon2] + 1.8)) ** exponent
#                 #if probCorrectionsPlus[phon1][phon2] > 1:
#                 #    correctionFactor *= 1 + log(probCorrectionsPlus[phon1][phon2]) ** exponent
#                 print phon1 + phon2 + ": " + str(probCorrectionsPlus[phon1][phon2]) + ", " + str(probCorrectionsMinus[phon1][phon2]) + " -> " + str(correctionFactor)
#     print
#     for phon1 in sounds:
#         for phon2 in sounds:
#             if probCorrectionsPlus[phon1][phon2] != 0 or probCorrectionsMinus[phon1][phon2] != 0:
#                 posEvidence = probCorrectionsPlus[phon1][phon2]
#                 if (posEvidence < 0.1): posEvidence = 0.1
#                 negEvidence = probCorrectionsMinus[phon1][phon2]
#                 if (negEvidence < 0.1): negEvidence = 0.1
#                 correctionFactor = (posEvidence / negEvidence) ** exponent
#                 #correctionFactor = 1.0
#                 #if probCorrectionsMinus[phon1][phon2] > 1:
#                 #    correctionFactor *= 1/(log(probCorrectionsMinus[phon1][phon2] + 1.8)) ** exponent
#                 #if probCorrectionsPlus[phon1][phon2] > 1:
#                 #    correctionFactor *= 1 + log(probCorrectionsPlus[phon1][phon2]) ** exponent
#                 global_replacement_probabilities[phon1][phon2] *= correctionFactor
#         normalizationFactor = sum(global_replacement_probabilities[phon1].values())
#         if normalizationFactor != 0.0:
#             for phon2 in sounds:
#                 global_replacement_probabilities[phon1][phon2] /= normalizationFactor
    
    replacementWeightsFile = open("replacement-weights" + str(iteration) + ".txt",'w')
    replacementWeightsFile.write("\t".join(sounds) + "\n")
    for phon1 in sounds:
        replacementWeightsFile.write("\t".join([str(- log(global_replacement_probabilities[phon1][phon2])) for phon2 in sounds]) + "\n")
    replacementWeightsFile.close()
    iteration += 1
