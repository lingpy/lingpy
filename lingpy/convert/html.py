# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-10-10 16:31
# modified : 2013-10-15 12:51
"""
Basic functions for HTML-plots.
"""

__author__="Johann-Mattis List"
__date__="2013-10-15"


import os
import colorsys
import codecs
import webbrowser

from ..settings import rcParams

import numpy as np


def colorRange(
        number,
        brightness = 300,
        ):
    """
    Function returns different colors for the given range.

    Notes
    -----
    Idea taken from
    http://stackoverflow.com/questions/876853/generating-color-ranges-in-python .
    """

    # get the hsv codes for the given range
    hsv = [(x * 1.0 / number, 0.5, 0.5) for x in range(number)]
    
    # convert hsv to rgb
    rgb = list(map(lambda x: colorsys.hsv_to_rgb(*x), hsv))
    
    # convert rgb to html
    for i in range(number):
        rgb[i] = tuple([int(rgb[i][k] * brightness) for k in range(3)])
    html = []
    for i in range(number):
        html.append("#%02x%02x%02x" % rgb[i])

    return html

def alm2html(
        infile,
        title = '',
        shorttitle = '',
        filename='',
        colored= False,
        verbose = True,
        show = True,
        main_template = '',
        table_template = '',
        dataset = '',
        **keywords
        ):
    """
    Convert files in ``alm``-format into colored ``html``-format.

    Parameters
    ----------

    title : str
        Define the title of the output file. If no title is provided, the
        default title ``LexStat - Automatic Cognate Judgments`` will be used.

    shorttitle : str
        Define the shorttitle of the ``html``-page. If no title is provided,
        the default title ``LexStat`` will be used.
    
    Notes
    -----
    The coloring of sound segments with respect to the sound class they belong
    to is based on the definitions given in the
    ``color`` :py:class:`~lingpy.data.model.Model`. It can easily be changed
    and adapted. 

    See also
    --------
    lingpy.output.plot.msa2html

    .. todo:: Format Change
       
       Adapt the input format in order to make it more flexible for general
       approaches, especially in collaboration with qlc.

    """
    # get the path to the templates
    #path = os.path.join(os.path.split(os.path.abspath(__file__))[0], '/templates/')
    
    # open the infile
    try:
        data = codecs.open(infile, "r", "utf-8").read()[:-1]
    except:
        data = codecs.open(infile+'.alm', "r", "utf-8").read()[:-1]

    # create the outfile
    if not filename:
        filename = rcParams['filename']
    
    # read in the templates
    path = os.path.dirname(os.path.realpath(__file__))
    
    if main_template:
        try:
            html = codecs.open(main_template,'r','utf-8').read()
        except:
            html = codecs.open(
                    os.path.join(path,'templates',main_template),
                    'r',
                    'utf-8'
                    ).read()
    else:
        html_path = os.path.join(path, 'templates', 'alm2html.html')
        html = codecs.open(html_path,'r','utf-8').read()
    
    if table_template:
        try:
            table = codecs.open(table_template,'r','utf-8').read()
        except:
            table = codecs.open(
                    os.path.join(path,'templates',table_template),
                    'r',
                    'utf-8'
                    ).read()
    else:
        table_path = os.path.join(path, 'templates', 'alm2html.table.html')
        table = codecs.open(table_path,'r','utf-8').read()

    # split the data into blocks
    blocks = data.split('\n\n')

    # retrieve the dataset
    if not dataset:
        dataset = blocks[0]

    # iterate over the rest of the blocks and store the data in a dictionary
    cogs = {}
    
    # create the outstring
    tmp_str = ''

    for block in blocks[1:]:
        lines = block.split('\n')
        
        m = []
        for l in lines:
            m.append(l.split('\t'))
        
        # create colordict for different colors
        dc = len(set([l[0] for l in m]))
        
        if colored:
            colors = dict([(a,b) for a,b in zip(
                sorted(set([int(l[0]) for l in m])),
                colorRange(
                    dc,
                    brightness = 400
                    ),
                )])
        else:
            colors = []
            white = True
            for i in sorted(set([abs(int(l[0])) for l in m])):
                if white:
                    colors.append((i,'white'))
                    white = False
                else:
                    colors.append((i,'lightgray'))
                    white = True
            colors = dict(colors)
        
        # get the basic item and its id
        iName = m[0][2]
        iID = m[0][3]

        # start writing the stuff to string
        tmp_str += table.format(
                NAME=iName,
                ID = iID)
        # define the basic string for the insertion
        bas = '<tr style="background-color:{0};font-weight:{2}">\n{1}\n</tr>'

        for tracer,l in enumerate(m):
            # check whether the current line is a borrowing
            if int(l[0]) < 0:
                loan_line = 'bold'
            else:
                loan_line = 'normal'

            # assign the cognate id
            tmp = '<td>{0}</td>\n'.format(l[0])
            tmp += '<td>{0}</td>\n'.format(l[1].strip('.'))
            tmp += '<td>{0}</td>\n'.format(''.join(l[4:]).replace('-',''))
            tmp += '<td style="background-color:{0}">'.format(colors[abs(int(l[0]))])
            tmp += '<table style="background-color:{0}">\n'.format(colors[abs(int(l[0]))])
            tmp += '<tr>\n{0}\n</tr>\n</table>\n'
            
            # check whether another entry follows that is also an alignment,
            # otherwise, there's no need to display a word as an alignment
            cognate_set = False
            if tracer < len(m)-1:
                if abs(int(m[tracer+1][0])) == abs(int(l[0])):
                    cognate_set = True
            if tracer > 0:
                if abs(int(m[tracer-1][0])) == abs(int(l[0])):
                    cognate_set = True

            if cognate_set: #len(l[4:]) > 1:
                alm = ''
                for char in l[4:]:
                    char = char
                    error = False
                    try:
                        c = rcParams['_color'][char]
                    except:
                        try:
                            c = rcParams['_color'][char[0]]
                        except:
                            c = 'white'
                            error = True
                    alm += '<td style="width:30px;text-align:center;'
                    if error:
                        alm += 'background-color:{0};color:red;font-weight:bold">{1}</td>'.format(c,char)
                       
                    else:
                        alm += 'background-color:{0};color:white;font-weight:bold;">{1}</td>'.format(c,char)
            else:
                alm = '<td style="border-color:{0};background-color:{1};">{0}'.format('--',colors[abs(int(l[0]))])
            
            # format the alignment
            tmp = tmp.format(alm)

            # check for last line, where a new line should be inserted (not the
            # fastest solution, but plotting is not a matter of time, and it
            # suffices it's current purpose
            if tracer < len(m)-1:
                pass
            else:
                tmp += '<tr style="background-color:white;"><td colspan="4">'
                tmp += '<hr style="background-color:white;border-color:white;border:0;height:2pt;"/>\n'
            
            # format the whole string
            tmp_str += bas.format(
                    colors[abs(int(l[0]))],
                    tmp,
                    loan_line
                    )

    if not title:
        title = "LexStat - Automatic Cognate Judgments"
    if not shorttitle:
        shorttitle = "LexStat"

    html = html.format(
            shorttitle = shorttitle,
            title = title,
            table = tmp_str,
            dataset = dataset,
            **keywords
            )

    out = codecs.open(filename+'.html','w','utf-8')
    out.write(html)
    out.close()

    if show:
        url = 'file://'+os.path.abspath(os.curdir)+'/'+filename+'.html'
        webbrowser.open(url)

    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.html'))
    return
