# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-09-15 21:41
# modified : 2013-09-15 21:41
"""
Module provides basic operations on Wordlist-Objects.
"""

__author__="Johann-Mattis List"
__date__="2013-09-15"

# external imports
import re
import json
import codecs

# internal imports
from ..settings import rcParams
from ..convert import *
from ..algorithm import clustering

def wl2dst(
        wl, # wordlist object
        taxa = "taxa",
        concepts = "concepts",
        ref = 'cogid',
        refB = '',
        mode = 'swadesh',
        **keywords
        ):
    """
    Function converts wordlist to distance matrix.
    """
    # check for taxon attribute
    taxa = getattr(wl,taxa)

    # check for concepts
    concepts = getattr(wl,concepts)

    distances = [[0 for i in range(wl.width)] for j in range(wl.width)]

    for i,taxA in enumerate(taxa):
        for j,taxB in enumerate(taxa):
            if i < j:
                if mode == 'shared':
                    listA = wl.get_list(col=taxA,entry=ref)
                    listB = wl.get_list(col=taxB,entry=ref)
                    
                    shared = [x for x in listA if x in listB]
                    missing = 0
                    score = len(shared)

                elif mode == 'swadesh':

                    # get the two dictionaries
                    dictA = wl.get_dict(col=taxA,entry=ref)
                    dictB = wl.get_dict(col=taxB,entry=ref)
    
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
                    try:
                        score = 1 - shared / (wl.height - missing)
                    except ZeroDivisionError:
                        print(rcParams['E_zero_division'].format(taxA,taxB))
                        score = 1.0
                
                 
                distances[i][j] = score
                if not refB:
                    distances[j][i] = score
            elif i == j:
                if mode == 'shared':
                    distances[i][j] = len(wl.get_list(col=taxA,flat=True))
            elif i > j and refB:
                if mode == 'shared':
                    listA = wl.get_list(col=taxA,entry=refB)
                    listB = wl.get_list(col=taxB,entry=refB)
                    
                    shared = [x for x in listA if x in listB]
                    missing = 0
                    score = len(shared)

                elif mode == 'swadesh':

                    # get the two dictionaries
                    dictA = wl.get_dict(col=taxA,entry=refB)
                    dictB = wl.get_dict(col=taxB,entry=refB)
    
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
                    try:
                        score = 1 - shared / (wl.height - missing)
                    except ZeroDivisionError:
                        print(rcParams['E_zero_division'].format(taxA,taxB))
                        score = 1.0
                 
                distances[i][j] = score

    return distances

def wl2dict(
        wordlist,
        sections,
        entries,
        exclude = []
        ):
    """
    Convert a wordlist to a complex dictionary with headings as keys.
    """
    
    # define output dictionary
    out = {}
    
    # determine the last section
    sorted_sections = sorted(sections)
    last_section = sorted_sections[-1]
   
    # iterate over wordlist
    for key in wordlist:
        
        if key not in exclude:
            # pass temporary pointer
            tmp = out

            # iterate over sections
            for s in sorted_sections:

                # get datapoint and text
                data_point = wordlist[key,sections[s][0]]
                text = sections[s][1].format(data_point)

                # dive deeper if this is not the last section
                if s != last_section:
                    
                    # access datapoint and text
                    #data_point = wordlist[key,sections[s][0]]
                    #text = sections[s][1].format(data_point)

                    # dive deeper
                    try:
                        tmp = tmp[data_point,text]
                    except KeyError:
                        tmp[data_point,text] = {}
                        tmp = tmp[data_point,text]
                
                # assign all entries
                else:
                    
                    # dive to last level
                    try:
                        tmp = tmp[data_point,text]
                    except:
                        tmp[data_point,text] = []
                        tmp = tmp[data_point,text]
                    
                    # get the final list of entries
                    tmp_list = []
                    for entry,format_string in entries:
                        if type(entry) in (list,tuple):
                            entry = ' '.join(entry)
                        tmp_list += [
                                format_string.format(
                                    wordlist[key,entry]
                                    )
                                ]
                    tmp += [tmp_list]
    
    return out

def renumber(wordlist,source,target=''):
    """
    Create numerical identifiers from string identifiers.
    """

    # iterate over wordlist and get all source ids
    sources = []
    for k in wordlist:
        sources += [wordlist[k,source]]

    sources = sorted(set(sources))

    # convert to numbers
    targets = list(range(1,len(sources)+1))

    # add to wordlist
    if not target:
        target = source+'id'

    # make converter
    converter = dict(zip(sources,targets))
    
    wordlist.add_entries(target,source,lambda x:converter[x])

    # add stuff to meta
    wordlist._meta[source+'2'+target] = converter

    if rcParams['verbose']: print("[i] Successfully renumbered {0}.".format(source))

def clean_taxnames(
        wordlist,
        taxa = 'taxa'
        ):
    # clean the names for all taxa in a wordlist
    current_taxa = eval('wordlist.'+taxa)
    new_taxa = [''.join([t for t in taxon if t not in ['() ,;[]{}']]) for taxon in current_taxa]

    old2new = dict(zip(current_taxa,new_taxa))

    wordlist.cols = [old2new[t] for t in wordlist.cols]
    
    # get idx for taxon
    idx = wordlist._header[taxa]

    for key in wordlist:
        wordlist._data[key][idx] = old2new[wordlist[key][idx]]
    wordlist.add_entries('taxa','taxa',lambda x:old2new[x],override=True)
    wordlist._clean_cache()
    return wordlist


def calculate(
        wordlist,
        data,
        taxa = 'taxa',
        concepts = 'concepts',
        ref = 'cogid',
        **keywords
        ):
    """
    Manipulate a wordlist object by adding different kinds of data.
    
    Parameters
    ----------
    data : str
        The type of data that shall be calculated. Currently supports

        * "tree": calculate a reference tree based on shared cognates
        * "dst": get distances between taxa based on shared cognates
        * "cluster": cluster the taxa into groups using different methods


    """
    defaults = dict(
            distances = False,
            tree_calc = "upgma",
            cluster = "upgma",
            force = False,
            threshold = 0.5,
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # get taxa for current calculation
    these_taxa = eval('wordlist.'+taxa)
    
    # calculate distances
    if data in ['distances','dst']:
        wordlist._meta['distances'] = wl2dst(
                wordlist,
                taxa,
                concepts,
                ref,
                **keywords
                )
    elif data in ['tre','tree','nwk']:
        if 'distances' not in wordlist._meta:
            distances = wl2dst(wordlist,taxa,concepts,ref,**keywords)
        else:
            distances = wordlist._meta['distances']
        if 'tree' in wordlist._meta and not keywords['force']:
            print(rcParams['W_force'].format('Reference tree'))
            return
        wordlist._meta['tree'] = clustering.matrix2tree(
                distances,
                these_taxa,
                keywords['tree_calc'],
                keywords['distances']
                )

    elif data in ['groups','cluster']:
        if 'distances' not in wordlist._meta:
            distances = wl2dst(wordlist,taxa,concepts,ref,**keywords)
        else:
            distances = wordlist._meta['distances']
        if 'groups' in wordlist._meta and not keywords['force']:
            print(rcParams['W_force'].format('Distance matrix'))
            return
        wordlist._meta['groups'] = clustering.matrix2groups(
                keywords['threshold'],
                distances,
                these_taxa
                )

    if rcParams['verbose']: print("[i] Successfully calculated {0}.".format(data))
        
def wl2qlc(
        header,
        data,
        filename = '',
        formatter = 'concept',
        **keywords
        ):
    """
    Write the basic data of a wordlist to file.
    """
    formatter = formatter.upper()

    defaults = dict(
            fileformat = 'qlc',
            ignore = []
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    if not filename:
        filename = rcParams['filename']

    # create output string
    out = '# Wordlist\n'

    # write meta to file
    if "meta" in keywords:
        meta = keywords["meta"]
    else:
        meta = {}
    
    kvpairs = {}
    jsonpairs = {}
    msapairs = {}
    trees = {}
    distances = ''
    taxa = ''
    scorer = ''

    for k,v in meta.items():
        # simple key-value-pairs
        if type(v) in [str,int] or k == "tree":
            kvpairs[k] = v
        elif k == 'msa':
            # go a level deeper, checking for keys
            for ref in v:
                if ref not in msapairs:
                    msapairs[ref] = {}
                for a,b in v[ref].items():
                    msapairs[ref][a] = b
        elif k == 'distances':
            distances = matrix2dst(v,meta['taxa'])
        elif k == 'taxa':
            taxa = '\n'.join(meta['taxa'])
        elif k == 'trees':
            trees = ''
            for key,value in v.items():
                trees += '<tre id="{0}">\n{1}\n</tre>\n'.format(
                        key,
                        value
                        )
        elif k == 'scorer' and k not in keywords['ignore']:
            for key,value in v.items():
                scorer += '<{2} id="{0}">\n{1}</{2}>\n\n'.format(
                        key,
                        scorer2str(value),
                        k
                        )
        else:
            # check whether serialization works
            try:
                json.dumps(v)
                jsonpairs[k] = v
            except TypeError:
                pass

    if kvpairs:
        out += '\n# META\n'
        for k,v in sorted(kvpairs.items(),key=lambda x:x[0]):
            out += '@{0}:{1}\n'.format(k,v)
    if taxa and keywords['taxa']:
        out += '\n# TAXA\n<taxa>\n'+taxa+'\n</taxa>\n'
    if jsonpairs:
        out += '\n# JSON\n'
        out += "<json>\n"
        out += json.dumps(jsonpairs,indent=4)
        out += '\n</json>\n'
    if msapairs:
        for ref in msapairs:
            out += "\n# MSA reference: {0}\n".format(ref)
            for k,v in msapairs[ref].items():
                
                if 'consensus' in v:
                    out += '#\n<msa id="{0}" ref="{1}" consensus="{2}">\n'.format(
                            k,ref,v['consensus'])
                else:
                    out += '#\n<msa id="{0}" ref="{1}">\n'.format(k,ref)
                out += msa2str(v,wordlist=True)
                out += "</msa>\n"
    if distances and 'distances' not in keywords['ignore']:
        out += '\n# DISTANCES\n<dst>\n'
        out += distances+'</dst>\n'

    if trees:
        out += '\n# TREES\n'+trees

    if scorer and 'scorer' not in keywords['ignore']:
        out += '\n# SCORER\n'+scorer

    out += '\n# DATA\nID\t'+'\t'.join(header)+'\n'
    
    # check for gloss in header to create nice output format
    if formatter in header:
        idx = header.index(formatter)
        formatter = None
    else:
        idx = False
        formatter = ''

    for key in sorted(data.keys()):
        
        # get the line
        line = data[key]
        
        # check for formatter
        if idx in range(len(line)):
            if line[idx] != formatter:
                out += '#\n'
                formatter = line[idx]

        # add the key 
        out += str(key)
        
        # add the rest of the values
        for value in line:
            if type(value) == list:
                try:
                    out += '\t'+' '.join(value)
                except:
                    out += '\t'+' '.join([str(v) for v in value])
            elif type(value) == int:
                out += '\t'+str(value)
            elif type(value) == float:
                out += '\t{0:.4f}'.format(value)
            else:
                out += '\t'+value
        out += '\n'

    f = codecs.open(filename +'.'+ keywords['fileformat'],'w','utf-8')
    f.write(out)
    if "stamp" in keywords:
        f.write(keywords['stamp'])
    f.close()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.'+keywords['fileformat']))

    return                 

def wl2csv(
        header,
        data,
        filename = '',
        formatter = 'concept',
        verbose = True,
        **keywords
        ):
    
    print("[WARNING] wl2csv is deprecated")
    return wl2qlc(header,data,filename,formatter,**keywords)
