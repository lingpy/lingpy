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

