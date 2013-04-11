# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-01-28 11:47
# modified : 2013-03-14 10:47
"""
Module provides functions for the transformation of text data into visually appealing format.
 
The main idea is to render alignments in colored tables where the colors of
the cells are chosen with respect to the sound-class of the sound value of each
given cell, such as in the following example: 

.. image:: colored_alignment.jpg
   :width: 600px
   
Here, the coloring of sounds follows the sound-class model of
:evobib:`Dolgopolsky1986`, the black margin around the cells in the second, the
third, and the fourth column indicates a swapped site.  The benefit of this way
to display alignments is that differences and similarities between the
sequences become visible at once, making it easy to check the correctness of a
given alignment analysis.
"""

__author__="Johann-Mattis List"
__date__="2013-03-14"

import os
import colorsys

from lingpy.check.exceptions import ThirdPartyModuleError

try:
    import networkx as nx
except ImportError:
    ThirdPartyModuleError('networkx').warning
 
from ..data import *
from ..data import _color
from ..align.sca import SCA
from ..check.messages import *


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
        title = None,
        shorttitle = None,
        filename=None,
        colored=True
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
    path = os.path.split(os.path.abspath(__file__))[0] + '/templates/'
    
    # open the infile
    try:
        data = open(infile).read()[:-1]
    except:
        data = open(infile+'.alm').read()[:-1]

    # create the outfile
    if not filename:
        outfile = infile.strip('.alm')+'.html'
    else:
        if filename.endswith('.html'):
            outfile = filename
        else:
            outfile = filename+'.html'
    
    # read in the templates
    html = open(path+'alm2html.html').read()
    table = open(path+'alm2html.table.html').read()

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
                        c = _color[char]
                    except:
                        try:
                            c = _color[char[0]]
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

    out = open(outfile,'w')
    out.write(html)
    out.close()

    FileWriteMessage(filename,'html').message('written')
    return


def patchy_alm2html(
        infile,
        title = None,
        shorttitle = None,
        filename=None,
        colored=True
        ):
    """
    Convert files in ``alm``-format into colored ``html``-format.
    """
    # get the path to the templates
    path = os.path.split(os.path.abspath(__file__))[0] + '/templates/'
    
    print("[i] {0}".format(infile))
    # open the infile
    try:
        data = open(infile).read()[:-1]
    except:
        data = open(infile+'.alm.patchy').read()[:-1]

    # create the outfile
    if not filename:
        outfile = infile.replace('.alm.patchy','')+'_patchy.html'
    else:
        if filename.endswith('_patchy.html'):
            outfile = filename
        else:
            outfile = filename+'_patchy.html'

    print("[i] "+outfile)

    
    # read in the templates
    html = open(path+'alm2html.html').read()
    table = open(path+'alm2html.table.html').read()

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
                sorted(set([float(l[0]) for l in m])),
                colorRange(
                    dc,
                    brightness = 400
                    ),
                )])
        else:
            colors = []
            white = True
            for i in sorted(set([float(l[0]) for l in m])):
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
                # this is the difference, it is not the last item that will be
                # accessed
                for char in l[4:]:
                    char = char
                    error = False
                    try:
                        c = _color[char]
                    except:
                        try:
                            c = _color[char[0]]
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
            # add the origin-key
            #tmp += '<td>{0}</td>'.format(l[-1])

            tmp += '<tr><td></td></tr>\n'
            tmp_str += bas.format(
                    colors[float(l[0])],
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

    out = open(outfile,'w')
    out.write(html)
    out.close()

def msa2tex(
        infile,
        template = None,
        path = None,
        filename = None
        ):
    """
    Convert an MSA to a tabular representation which can easily be used in
    LaTeX documents.
    """
    # while alm-format can be read from the text-file without problems,
    # msa-format should be loaded first (once this is already provided), the
    # loss in speed won't matter much, since output of data is not a daily task
    
    ## get the path to the templates
    if not path:
        path = os.path.split(os.path.abspath(__file__))[0] + '/templates/'
    else:
        if path.endswith('/'):
            pass
        else:
            path += '/'

    # load msa
    msa = SCA(infile)

    ## load templates
    if not template:
        tex = open(path+'msa.tex').read()
    else:
        tex = open(template).read()

    # load dataset, etc.

    dataset = msa.dataset
    pid_score = int(100 * msa.get_pid(2) + 0.5)
    infile = msa.infile
    seq_id = msa.seq_id
    
    # determine the length of the longest taxon
    taxl = max([len(t) for t in msa.taxa])
    
    height = len(msa.alm_matrix)
    width = len(msa.alm_matrix[0])
    
    start = r'\tabular{l'+width*'c'+'}\n'
    start += r'\bf\ttfamily Taxon & \multicolumn{'+str(width)+r'}{l}{\bf\ttfamily Alignment}\\'+'\n'

    # get the dolgo model for colors
    msa.ipa2cls(dolgo)
    
    # check for swaps in the alignment
    if hasattr(msa,'swap_index'):
        swaps = []
        for s in msa.swap_index:
            swaps.extend(s)
    else:
        swaps = []

    body = start
    for i,taxon in enumerate(msa.taxa):
        body += r'\ttfamily '+taxon.replace('_',r'\_')
        for j,cls in enumerate(msa.classes[i]):
            char = msa.alm_matrix[i][j]
            if char == '_':
                char = r'\#'
            if cls == '_':
                cls = '2'
            if j not in swaps:
                body += r'&\cellcolor{col'+cls+r'}'+char
            else:
                if char != '-':
                    body += r'&\cellcolor{col'+cls+r'}\color{white}\bf '+char
                else:
                    body += r'&\cellcolor{col'+cls+r'}\bf '+char
        body += r'\\'+'\n'

    body += r'&'+'&'.join([r'\color{white}XXX' for i in range(width)])+r'\\'+'\n'
    body += r'\endtabular'+'\n'
    
    # create the parameters etc.
    w = 1.5 * width + taxl * 0.25
    h = 0.5 * height + 1.0

    tex = tex.replace('<+WIDTH+>','{0:2f}'.format(w))
    tex = tex.replace('<+HEIGHT+>','{0:2f}'.format(h))

    # create the rput stuff
    tex = tex.replace('<+NEWX+>','{0:.2f}'.format(w/2.0))
    tex = tex.replace('<+NEWY+>','{0:.2f}'.format((h-0.5)/2.0))

    # insert the rest
    tex = tex.replace('<+CONTENT+>',body)

    # write to file
    if not filename:
        outfile = msa.infile + '.tex'
    else:
        outfile = filename + '.tex'

    out = open(outfile,'w')
    out.write(tex)
    out.close()


def string2html(
        taxon,
        string,
        swaps = [],
        tax_len = None,
        path = None,
        template = None
        ):
    """
    Function converts an (aligned) string into colored html-format.

    """

    # determine the length of the string
    if not tax_len:
        tax_len = len(taxon)

    # set the tr-line
    tr = '<tr class="msa">\n{0}\n</tr>'

    # set the td_taxon-line
    td_taxon = '<td class="taxon" width="'+str(15 * tax_len)+'">{0}</td>\n'

    # get the percentage scaling factor
    perc = int(80 / len(string) + 0.5)

    # get vals for residue and swaps
    td_residue = '<td class="residue" width="50" align="center" bgcolor="{1}">'+\
            '<font color="{2}">{0}</font></td>\n'
    td_swap = '<td class="residue swap" style="border:solid 3px black" width="50"'+\
            'align="center" bgcolor="{1}"><font color="{2}">{0}</font></td>\n'

    # start with filling the taxon
    out = ''
    out += td_taxon.format(taxon)
    
    # go on with the colors
    for i,char in enumerate(string):
        try:
            c = _color[char]
            fg = '#000000'
        except:
            try:
                c = _color[char[0]]
                fg = '#000000'
            except:
                input("Unknown character '"+char+"', press ANY key to continue. " )
                c = '#ffffff'
                fg = '#eb3410'

        if i in swaps:
            out += td_swap.format(char,c,fg)
        else:
            out += td_residue.format(char,c,fg)

    return out

def msa2html(
        infile,
        shorttitle = None,
        filename = None,
        path = None,
        template = None
        ):
    """
    Convert files in ``msa``-format into colored ``html``-format.

    Parameters
    ----------

    shorttitle : str
        Define the shorttitle of the ``html``-page. If no title is provided,
        the default title ``SCA`` will be used.

    filename : str
        Define the name of the output file. If no name is defined, the name of
        the input file will be taken as a default.

    Examples
    --------
    Load the libary.

    >>> from lingpy import *
    
    Load an ``msq``-file from the test-sets.

    >>> msa = Multiple(get_file('test.msq'))

    Align the data progressively and carry out a check for swapped sites.

    >>> msa.prog_align()
    >>> msa.swap_check()
    >>> print(msa)
    w    o    l    -    d    e    m    o    r    t
    w    a    l    -    d    e    m    a    r    -
    v    -    l    a    d    i    m    i    r    -

    Save the data to the file ``test.msa``.

    >>> msa.output('msa')

    Convert the ``msa``-file to ``html``.

    >>> msa2html('test.msa')
    
    Notes
    -----
    The coloring of sound segments with respect to the sound class they belong
    to is based on the definitions given in the
    ``color`` :py:class:`~lingpy.data.model.Model`. It can easily be changed
    and adapted. 

    See also
    --------
    lingpy.output.plot.alm2html
    """
    
    # while alm-format can be read from the text-file without problems,
    # msa-format should be loaded first (once this is already provided), the
    # loss in speed won't matter much, since output of data is not a daily task
    
    # get the path to the templates
    if not path:
        path = os.path.split(os.path.abspath(__file__))[0] + '/templates/'
    else:
        if path.endswith('/'):
            pass
        else:
            path += '/'

    # load msa
    msa = SCA(infile)

    # load templates
    if not template:
        html = open(path+'msa2html.html').read()
    else:
        html = open(template).read()

    # load dataset, etc.
    dataset = msa.dataset
    pid_score = int(100 * msa.get_pid(2) + 0.5)
    infile = msa.infile
    seq_id = msa.seq_id

    # define the titles etc.
    if not shorttitle:
        shorttitle = 'SCA'
    
    # determine the length of the longest taxon
    taxl = max([len(t) for t in msa.taxa])

    out = ''
    tr = '<tr class="msa">\n{0}\n</tr>'
    td_taxon = '<td class="taxon" width="'+str(15 * taxl)+'">{0}</td>\n'
    perc = int(80 / len(msa.alm_matrix[0]) + 0.5)
    td_residue = '<td class="residue" width="50" align="center" bgcolor="{1}">'+\
            '{0}</td>\n'
    td_swap = '<td class="residue swap" style="border:solid 3px black" width="50"'+\
            'align="center" bgcolor="{1}">{0}</td>\n'
    
    # check for swaps in the alignment
    if hasattr(msa,'swap_index'):
        swaps = []
        for s in msa.swap_index:
            swaps.extend(s)
    else:
        swaps = []

    # start iteration
    for i,taxon in enumerate(msa.taxa):
        tmp = ''
        tmp += td_taxon.format(taxon)
        for j,char in enumerate(msa.alm_matrix[i]):
            try:
                c = _color[char]
            except:
                try:
                    c = _color[char[0]]
                except:
                    print(char)
                    #input()            
            if j in swaps:
                tmp += td_swap.format(char,c)
            else:
                tmp += td_residue.format(char,c)
        out += tr.format(tmp)
    html = html.format(
            table = out,
            dataset = dataset,
            pid = pid_score,
            file = infile,
            sequence = seq_id,
            shorttitle = shorttitle,
            width=len(msa.alm_matrix[0]),
            table_width='{0}'.format(len(msa.alm_matrix[0])* 50 + 15 * taxl),
            taxa = len(msa.alm_matrix),
            uniseqs=len(msa.uniseqs)
            )
    
    
    if not filename:
        outfile = msa.infile + '.html'
    else:
        outfile = filename + '.html'

    # check, whether the outfile already exists
    out = open(outfile,'w')
    out.write(html)
    out.close()
