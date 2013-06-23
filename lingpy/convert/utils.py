# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-06-11 21:02
# modified : 2013-06-11 21:02
"""
Module provides various functions for specific conversion tasks.
"""

__author__="Johann-Mattis List"
__date__="2013-06-11"


def wl2dict(
        wordlist,
        sections,
        entries
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


#def wl2file(
#        wordlist,
#        sections = {},
#        entries = [],
#        filename = 'output',
#        fileformat = 'txt',
#        entry_sep = '',
#        #section_sep = '\n',
#        item_sep = '',
#        template = ''
#        ):
#    """
#    Export a wordlist to file with specific headings.
#    """
#    # check for sections
#    if not sections:
#        sections = dict(
#                h1 = ('concept','\n# Concept: {0}\n'),
#                h2 = ('cogid','## Cognate-ID: {0}\n'),
#                )
#
#    # check for entries
#    if not entries:
#        entries = [
#                ('language','{0}'),
#                ('counterpart',' {0} '),
#                ('ipa','{0}\n'),
#
#                ]
#    
#    # get the temporary dictionary
#    out = wl2dict(
#            wordlist,
#            sections,
#            entries
#            )
#
#    # assign the output string
#    out_string = ''
#
#    # iterate over the dictionary and start to fill the string
#    for key in sorted(out):
#        
#        # write key to file
#        out_string += key[1]
#        
#        # reassign tmp
#        tmp = out[key]
#
#        # set the pointer and the index
#        pointer = []
#        idx = 0
#
#        while True:
#            
#            # dive deeper
#            break_loop = False
#            
#            # check for type of current point
#            if type(tmp) == dict:
#                
#                # set the pointer first val points to the current depth, second
#                # to the keys which are successively popped
#                pointer += [[tmp,sorted(tmp.keys())]]
#                try:
#                    next_key = pointer[idx][1].pop(0)
#                    out_string += next_key[1]
#                    tmp = pointer[idx][0][next_key] 
#                    idx += 1
#                except: 
#                    break_loop = True
#            else:
#                tmp_strings = []
#                for line in sorted(tmp):
#                    tmp_strings += [item_sep.join(line)]
#                out_string += entry_sep.join(tmp_strings)
#                idx -= 1
#                tmp = pointer[idx][0]
#
#            if break_loop:
#                break
#    
#    # load the template
#    if template:
#        tmpl = open(template,'r').read()
#    else:
#        tmpl = '{0}'
#    
#    # open outfile
#    f = open(filename+'.'+fileformat,'w')
#    if fileformat == 'tex':
#        f.write(tmpl.format(out_string.replace('_',r'\_')))
#    else:
#        f.write(tmpl.format(out_string))
#    f.close()

