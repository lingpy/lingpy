"""
Module provides basic operations on Wordlist-Objects.
"""
from __future__ import unicode_literals, print_function, absolute_import, division

# external imports
import json
import unicodedata

from six import text_type

# internal imports
from ..settings import rcParams
from ..convert.strings import matrix2dst, scorer2str, msa2str, pap2nex, pap2csv
from ..algorithm import clustering
from .. import util
from .. import log


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
    logger = log.get_logger()
    # check for taxon attribute
    taxa = getattr(wl,taxa)

    # check for concepts
    concepts = getattr(wl,concepts)

    distances = [[0 for i in range(wl.width)] for j in range(wl.width)]

    for i,taxA in enumerate(taxa):
        for j,taxB in enumerate(taxa):
            if i < j:
                if mode in ['shared', 'jaccard']:
                    listA = wl.get_list(col=taxA,entry=ref)
                    listB = wl.get_list(col=taxB,entry=ref)
                    
                    shared = [x for x in listA if x in listB]
                    missing = 0

                    if mode == 'jaccard':
                        allc = len(set(listA+listB))
                        score = 1 - len(set(shared)) / allc
                    else:
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
                        logger.exception(
                            "Zero-division error encountered in '{0}' and '{1}'.".format(
                               taxA, taxB))
                        score = 1.0
                
                 
                distances[i][j] = score
                if not refB:
                    distances[j][i] = score
            elif i == j:
                if mode == 'shared':
                    distances[i][j] = len(wl.get_list(col=taxA,flat=True))
            elif i > j and refB:
                if mode in ['shared', 'jaccard']:
                    listA = wl.get_list(col=taxA,entry=refB)
                    listB = wl.get_list(col=taxB,entry=refB)
                    
                    shared = [x for x in listA if x in listB]
                    missing = 0

                    if mode == 'jaccard':
                        allc = len(set(listA+listB))
                        score = 1 - len(set(shared)) / allc
                    else:
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
                        logger.exception(
                            "Zero-division error encountered in '{0}' and '{1}'.".format(
                                taxA, taxB))
                        score = 1.0
                 
                distances[i][j] = score

    return distances

def wl2dict(
        wordlist,
        sections,
        entries,
        exclude = None
        ):
    """
    Convert a wordlist to a complex dictionary with headings as keys.
    """
    
    # define output dictionary
    out = {}
    
    if not exclude:
        exclude = []

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

def renumber(wordlist,source,target='', override=False):
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
    
    wordlist.add_entries(target,source,lambda x:converter[x], override=override)

    # add stuff to meta
    wordlist._meta[source+'2'+target] = converter

    log.info("Successfully renumbered {0}.".format(source))


def clean_taxnames(
        wordlist,
        column = 'doculect',
        f = lambda x:''.join([t for t in x if t not in '()[]{},;:'])
        ):
    """
    Function cleans taxon names in order to make sure they can be used in Newick files.

    """
    # clean the names for all taxa in a wordlist
    current_taxa = eval('wordlist.'+column)
    new_taxa = [f(taxon) for taxon in current_taxa]

    old2new = dict(zip(current_taxa,new_taxa))
    new2old = dict(zip(new_taxa,current_taxa))
    
    if column == wordlist._col_name:
        wordlist.cols = [old2new[t] for t in current_taxa]
    
    wordlist.add_entries('_doculect','doculect',lambda x:old2new[x],override=True)
    wordlist.add_entries('doculect','_doculect',lambda x:new2old[x],override=True)

def calculate_data(
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
    logger = log.get_logger()
    defaults = dict(
            distances = False,
            tree_calc = "upgma",
            cluster = "upgma",
            force = False,
            threshold = 0.5,
            cluster_method = 'upgma',
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
            wordlist._meta['distances'] = wl2dst(wordlist,taxa,concepts,ref,**keywords)
            distances = wordlist._meta['distances']
        else:
            distances = wordlist._meta['distances']
        if 'tree' in wordlist._meta and not keywords['force']:
            logger.warn("Reference tree has already been calculated, force overwrite by setting 'force' to 'True'.")
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
            logger.warn("Distance matrix has already been calculated, force overwrite by setting 'force' to 'True'.")
            return
        wordlist._meta['groups'] = clustering.matrix2groups(
                keywords['threshold'],
                distances,
                these_taxa,
                keywords['cluster_method']
                )
    log.info("Successfully calculated {0}.".format(data))


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
    defaults = dict(
            ignore = ['taxa', 'doculects', 'msa']
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    if keywords['ignore'] == 'all':
        keywords['ignore'] = ['taxa', 'doculects', 'msa', 'json']

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
        if type(v) in [text_type,int] or k == "tree":
            kvpairs[k] = v
        elif k == 'msa' and k not in keywords['ignore']:
            # go a level deeper, checking for keys
            for ref in v:
                if ref not in msapairs:
                    msapairs[ref] = {}
                for a,b in v[ref].items():
                    msapairs[ref][a] = b
        elif k == 'distances':
            distances = matrix2dst(v,meta['taxa'])
        elif k in ['taxa', 'doculect', 'taxon', 'doculects']: # and k not in keywords['ignore']:
            # we need to find a better solution here, since it is not nice to
            # have taxa written to json again and again
            pass
            #taxa = '\n'.join(meta['taxa'])
        elif k == 'trees' and k not in keywords['ignore']:
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
    if jsonpairs and 'json' not in keywords['ignore']:
        out += "@json: "+json.dumps(jsonpairs)+'\n'
        #out += '\n# JSON\n'
        #out += "<json>\n"
        #out += json.dumps(jsonpairs,indent=4)
        #out += '\n</json>\n'
    if msapairs and 'msa' not in keywords['ignore']:
        for ref in msapairs:
            out += "\n# MSA reference: {0}\n".format(ref)
            for k,v in msapairs[ref].items():
                if 'consensus' in v:
                    out += '#\n<msa id="{0}" ref="{1}" consensus="{2}">\n'.format(
                            k,ref,' '.join(v['consensus']))
                else:
                    out += '#\n<msa id="{0}" ref="{1}">\n'.format(k,ref)
                outs = msa2str(v,wordlist=True)
                out += outs
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
        sorted_data = sorted(data.keys(),key=lambda x:data[x][idx])
    elif len(formatter.split(',')) == 2:
        idxA,idxB = formatter.split(',')
        idxA = header.index(idxA)
        idxB = header.index(idxB)
        idx = idxA
        formatter = None
        sorted_data = sorted(data.keys(),key=lambda x:(data[x][idxA],data[x][idxB]))
    else:
        idx = False
        formatter = ''
        sorted_data = sorted(data.keys())

    for key in sorted_data:
        
        # get the line
        line = data[key]
        
        # check for formatter
        if idx in range(len(line)):
            if line[idx] != formatter:
                out += '#\n'
                formatter = line[idx]

        # add the key 
        out += text_type(key)
        
        # add the rest of the values
        for value in line:
            if type(value) == list:
                try:
                    out += '\t'+' '.join(value)
                except:
                    out += '\t'+' '.join([text_type(v) for v in value])
            elif type(value) == int:
                out += '\t'+text_type(value)
            elif type(value) == float:
                out += '\t{0:.4f}'.format(value)
            else:
                out += '\t'+value
        out += '\n'

    path = filename + '.' + keywords['fileformat']
    util.write_text_file(
        path, unicodedata.normalize("NFC", out) + keywords.get('stamp', ''))
    return

def wl2csv(
        header,
        data,
        filename = '',
        formatter = 'concept',
        verbose = True,
        **keywords
        ):
    log.deprecated('wl2csv', '')
    return wl2qlc(header,data,filename,formatter,**keywords)


def tsv2triple(wordlist, outfile=None):
    """
    Function converts a wordlist to a triple data structure.

    Notes
    -----
    The basic values of which the triples consist are:
      * ID (the ID in the TSV file)
      * COLUMN (the column in the TSV file)
      * VALUE (the entry in the TSV file)
    """
    
    tstore = []
    for head in wordlist.header:
        log.debug('tsv2triple: ' + head)
        for key in wordlist:
            tstore += [(key,head.upper(),wordlist[key,head])]

    if outfile:
        out = ''
        for a, b, c in tstore:
            if type(c) == list:
                c = ' '.join([text_type(x) for x in c])
            if c != '-':
                out += '{0}\t{1}\t{2}\n'.format(a, b, c)
        util.write_text_file(outfile, out, normalize='NFC')
    else:
        return tstore

def triple2tsv(infile, output="table"):
    """
    Function reads in a triple file and converts it to a tabular data structure.
    """
    
    D = {}
    idxs = set([])
    cols = set([])
    
    if not isinstance(infile, list):
        lines = util.read_text_file(infile, normalize='NFD', lines=True)
    else:
        lines = infile # shallow copy is ok here


    for line in lines: 
        a,b,c = line.split('\t')
        try:
            D[a][b] = c
        except KeyError:
            try:
                D[a] = {b:c}
            except KeyError:
                D = {a:{b:c}}

        idxs.add(a)
        cols.add(b)
    
    idxs = sorted(idxs)
    cols = sorted(cols)
    table = [[0] + ['' for col in cols] for idx in idxs]
    
    for i,idx in enumerate(idxs):
        table[i][0] = idx
        for j,col in enumerate(cols):
            try:
                table[i][j+1] = D[idx][col]
            except:
                pass
    
    output_table = [["ID"] + cols]+table

    if output in ['wordlist', 'dict']:
        wlD = {}
        for line in table:
            wlD[int(line[0])] = line[1:]
        wlD[0] = cols 

        return wlD
    
    else:
        return output_table
