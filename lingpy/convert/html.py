# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-10-10 16:31
# modified : 2013-10-10 16:31
"""
Basic functions for HTML-plots.
"""

__author__="Johann-Mattis List"
__date__="2013-10-10"


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
        colored=True,
        verbose = True,
        show = True
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
    html_path = os.path.join(path, 'templates', 'alm2html.html')
    table_path = os.path.join(path, 'templates', 'alm2html.table.html')
    html = codecs.open(html_path,'r','utf-8').read()
    table = codecs.open(table_path,'r','utf-8').read()

    # split the data into blocks
    blocks = data.split('\n\n')

    # retrieve the dataset
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
            for i in sorted(set([int(l[0]) for l in m])):
                if white:
                    colors.append((i,'white'))
                    white = False
                else:
                    colors.append((i,'gray'))
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
        bas = '<tr bgcolor="{0}">\n{1}\n</tr>'
        for l in m:
            # assign the cognate id
            tmp = '<td>{0}</td>\n'.format(l[0])
            tmp += '<td bgcolor=white> </td>\n'
            tmp += '<td>{0}</td>\n'.format(l[1].strip('.'))
            tmp += '<td bgcolor=white> </td>\n'
            tmp += '<td>{0}</td>\n'.format(''.join(l[4:]).replace('-',''))
            tmp += '<td bgcolor=white> </td>\n'
            tmp += '<td bgcolor=white><table bgcolor=white>\n<tr>\n{0}\n</tr>\n</table>\n'
            if len(l[4:]) > 1:
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
                    alm += '<td width="30px" align="center"'
                    if error:
                        alm += 'bgcolor="{0}"><font color="red"><b>{1}</b></font></td>'.format(c,char)
                       
                    else:
                        alm += 'bgcolor="{0}"><font color="white"><b>{1}</b></font></td>'.format(c,char)
            else:
                alm = '<td bgcolor="white">{0}'.format('--')

            tmp = tmp.format(alm)
            tmp += '<tr><td></td></tr>\n'
            tmp_str += bas.format(
                    colors[int(l[0])],
                    tmp)
    if not title:
        title = "LexStat - Automatic Cognate Judgments"
    if not shorttitle:
        shorttitle = "LexStat"

    html = html.format(
            shorttitle = shorttitle,
            title = title,
            table = tmp_str,
            dataset = dataset
            )

    out = codecs.open(filename+'.html','w','utf-8')
    out.write(html)
    out.close()

    if show:
        url = 'file://'+os.path.abspath(os.curdir)+'/'+filename+'.html'
        webbrowser.open(url)

    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.html'))
    return
