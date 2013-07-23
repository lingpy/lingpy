# author   : Johann-Mattis List
# email    : mattis.list@gmail.com
# created  : 2013-01-28 11:47
# modified : 2013-07-10 12:37
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
__date__="2013-07-10"

import os
import colorsys
import codecs

from ..settings import rcParams

try:
    import networkx as nx
except ImportError:
    print(rcParams['W_missing_module'].format('networkx'))
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except:
    print(rcParams['W_missing_module'].format('matplotlib'))

from ..align.sca import SCA
from ..thirdparty import cogent as cg
from .gml import *

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
        verbose = True
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
        filename = 'lingpy-{0}'.format(_timestamp())
    
    # read in the templates
    html = codecs.open(path+'alm2html.html','r','utf-8').read()
    table = codecs.open(path+'alm2html.table.html','r','utf-8').read()

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

    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.html'))
    return

def patchy_alm2html(
        infile,
        title = None,
        shorttitle = None,
        filename=None,
        colored=True,
        verbose = False
        ):
    """
    Convert files in ``alm``-format into colored ``html``-format.
    """
    # get the path to the templates
    path = os.path.split(os.path.abspath(__file__))[0] + '/templates/'
    
    if rcParams['verbose']: print("[i] {0}".format(infile))
    # open the infile
    try:
        data = open(infile).read()[:-1]
    except:
        data = open(infile+'.alm.patchy').read()[:-1]

    # create the outfile
    if not filename:
        filename = 'lingpy-{0}_patchy'.format(_timestamp())

    if rcParams['verbose']: print("[i] "+filename)
    
    # read in the templates
    html = codecs.open(path+'alm2html.html','r','utf-8').read()
    table = codecs.open(path+'alm2html.table.html','utf-8').read()

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

    out = codecs.open(filename+'.patchy.html','w','utf-8')
    out.write(html)
    out.close()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.patchy.html'))

def msa2tex(
        infile,
        template = '',
        path = '',
        filename = '',
        verbose = True
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
        tex = codecs.open(path+'msa.tex','r','utf-8').read()
    else:
        tex = codecs.open(template,'r','utf-8').read()

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
        filename = 'lingpy-{0}'

    out = codecs.open(filename+'.tex','w','utf-8')
    out.write(tex)
    out.close()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.tex'))


def string2html(
        taxon,
        string,
        swaps = [],
        tax_len = None,
        path = '',
        template = ''
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
            c = rcParams['_color'][char]
            fg = '#000000'
        except:
            try:
                c = rcParams['_color'][char[0]]
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
        shorttitle = '',
        filename = '',
        path = '',
        template = '',
        verbose = True
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
        html = codecs.open(path+'msa2html.html','r','utf-8').read()
    else:
        html = codecs.open(template,'r','utf-8').read()

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
                c = rcParams['_color'][char]
            except:
                try:
                    c = rcParams['_color'][char[0]]
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
        filename = rcParams['filename']

    # check, whether the outfile already exists
    out = codecs.open(outfile,'w','utf-8')
    out.write(html)
    out.close()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.html'))

def plot_gls(
        gls,
        treestring,
        degree = 90,
        fileformat = 'pdf',
        verbose = True,
        **keywords
        ):
    """
    Plot a gain-loss scenario for a given reference tree.
    """

    # get kewyords
    defaults = dict(
            figsize = (15,15),
            left = 0.05,
            top = 0.95,
            bottom = 0.05,
            right = 0.95,
            radius = 0.5,
            textsize = 8,
            edgewidth = 5,
            linewidth = 2,
            scale_radius = 1.2,
            ylim = 1,
            xlim = 1,
            text = True,
            gain_color = 'white',
            loss_color = 'black',
            gain_linestyle = 'dotted',
            loss_linestyle = 'solid',
            ax_linewidth = 0,
            filename = 'lingpy-{0}'.format(_timestamp)
            )

    for k in defaults: 
        if k not in keywords:
            keywords[k] = defaults[k]

    # set filename as variabel for convenience
    filename = keywords['filename']
    
    try:
        tree = cg.LoadTree(treestring=treestring)
    except:
        try:
            tree = cg.LoadTree(treestring)
        except:
            tree = treestring

    tgraph = radial_layout(treestring,degree=degree)
    
    graph = gls2gml(
            gls,
            tgraph,
            tree
            )

    nodes = []

    # assign nodes and edges
    for n,d in graph.nodes(data=True):
        g = d['graphics']
        x = g['x']
        y = g['y']
        s = d['state']
    
        nodes += [(x,y,s)]

    # now plot the stuff
    fig = plt.figure(figsize = keywords['figsize'])
    figsp = fig.add_subplot(111)
    figsp.axes.get_xaxis().set_visible(False)
    figsp.axes.get_yaxis().set_visible(False)

    # set the axes linewidht
    for s in figsp.spines.values():
        s.set_linewidth(keywords['ax_linewidth'])
    
    #plt.axes(frameon=keywords['frameon'])
    plt.axis('equal')

    for nA,nB in graph.edges():
        xA = graph.node[nA]['graphics']['x']
        xB = graph.node[nB]['graphics']['x']
        yA = graph.node[nA]['graphics']['y']
        yB = graph.node[nB]['graphics']['y']

        plt.plot(
                [xA,xB],
                [yA,yB],
                '-',
                color = 'black',
                linewidth = keywords['edgewidth'],
                zorder = 1
                )

    # now, iterate over nodes
    for x,y,s in nodes:
        if s == 'O':
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius'],
                    0,360,
                    facecolor = keywords['gain_color'],
                    linewidth = keywords['linewidth'],
                    linestyle = keywords['gain_linestyle']
                    )
        elif s == 'o':
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius'] / keywords['scale_radius'],
                    0,360,
                    facecolor = keywords['gain_color'],
                    linewidth = keywords['linewidth']
                    )
        elif s == 'L':
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius'],
                    0,360,
                    facecolor = keywords['loss_color'],
                    linewidth = keywords['linewidth'],
                    linestyle = keywords['loss_linestyle']
                    )
        else:
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius'] / keywords['scale_radius'],
                    0,360,
                    facecolor = keywords['loss_color'],
                    linewidth = keywords['linewidth']
                    )
        figsp.add_artist(w)

        # if text is chosen as argument
        if keywords['text']:
            if s in 'Oo': 
                t = '1'
                c = 'black'
            else: 
                t = '0'
                c = 'white'

            plt.text(
                    x,
                    y,
                    t,
                    size = keywords['textsize'],
                    color = c,
                    va = "center",
                    ha = "center",
                    fontweight = 'bold'
                    )

    # set x and y-values
    xvals = [x[0] for x in nodes]
    yvals = [x[1] for x in nodes]
    
    plt.xlim(min(xvals)-keywords['xlim'],max(xvals)+keywords['xlim'])
    plt.ylim(min(yvals)-keywords['ylim'],max(yvals)+keywords['ylim'])
    
    plt.subplots_adjust(
            left = keywords['left'],
            right = keywords['right'],
            top = keywords['top'],
            bottom = keywords['bottom']
            )
    plt.savefig(
            filename + '.'+fileformat
            )
    plt.clf()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.'+fileformat))

def plot_tree(
        treestring,
        degree = 90,
        fileformat = 'pdf',
        verbose = True,
        **keywords
        ):
    """
    Plot a newick tree to PDF or other graphical formats.
    """

    default = dict(
            linewidth = 5,
            linecolor = 'black',
            nodesize = 10,
            nodecolor = 'black',
            textsize = '10',
            textcolor = 'white',
            va = 'center',
            ha = 'center',
            bg = 'black',
            fontweight = 'bold',
            left = 0.05,
            right = 0.95,
            top = 0.95,
            bottom = 0.05,
            figsize = (10,10),
            node_dict = {},
            no_labels = False,
            xlim            = 5,
            ylim            = 5,
            xlimr           = False,
            xliml           = False,
            ylimt           = False,
            ylimb           = False,
            change = lambda x: x**1.75,
            frameon = False,
            edge_list = [],
            ax_linewidth = 0,
            start = 0,
            usetex = False,
            filename = 'lingpy-{0}'.format(_timestamp())
            )
    for k in default:
        if k not in keywords:
            keywords[k] = default[k]

    # set filename as variable for convenience
    filename = keywords['filename']

    # switch backend, depending on whether tex is used or not
    backend = mpl.get_backend()
    if keywords['usetex'] and backend != 'pgf':
        plt.switch_backend('pgf')
    elif not keywords['usetex'] and backend != 'TkAgg':
        plt.switch_backend('TkAgg')

    
    # get the tree-graph
    graph = radial_layout(
            treestring,
            degree=degree,
            change=keywords['change'],
            start=keywords['start']
            )
    
    # create the figure
    fig = plt.figure(figsize=keywords['figsize'])
    figsp = fig.add_subplot(111)
    figsp.axes.get_xaxis().set_visible(False)
    figsp.axes.get_yaxis().set_visible(False)

    for s in figsp.spines.values():
        s.set_linewidth(keywords['ax_linewidth'])

    #plt.axes(frameon=keywords['frameon'])
    plt.axis('equal')
    plt.xticks([])
    plt.yticks([])
    
    # get xlim and ylim
    xvals,yvals = [],[]
    # start iterating over edges
    for nA,nB,d in graph.edges(data=True)+keywords['edge_list']:
            
        # get the coordinates
        xA = graph.node[nA]['graphics']['x']
        yA = graph.node[nA]['graphics']['y']
        xB = graph.node[nB]['graphics']['x']
        yB = graph.node[nB]['graphics']['y']
        
        if 'color' in d:
            plt.plot(
                    [xA,xB],
                    [yA,yB],
                    '-',
                    **d
                    )
        else:
            plt.plot(
                    [xA,xB],
                    [yA,yB],
                    '-',
                    color = keywords['linecolor'],
                    linewidth = keywords['linewidth'],
                    )

    # get the nodes
    for n,d in graph.nodes(data=True):
           
        g = d['graphics']
        x,y = g['x'],g['y']

        xvals += [x]
        yvals += [y]

        # try to get information from the node-dict
        try:
            settings = keywords['node_dict'][n]
        except:
            settings = {}
        
        # overwrite the stuff in keywords
        for k in keywords:
            if k not in settings:
                settings[k] = keywords[k]
    
        if d['label'].startswith('edge') or d['label'].startswith('root') or keywords['no_labels']:
            plt.plot(
                    x,
                    y,
                    'o',
                    markersize = settings['nodesize'],
                    color = settings['nodecolor'],
                    markeredgewidth = settings['linewidth']
                    )
        else:
            plt.text(
                    x,
                    y,
                    d['label'],
                    color = settings['textcolor'],
                    fontweight = settings['fontweight'],
                    va = settings['va'],
                    ha = g['s'],
                    bbox = dict(
                        facecolor=settings['bg'],
                        boxstyle = 'square,pad=0.2',
                        ec="none",
                        ),
                    size = settings['textsize'],
                    rotation = g['angle'],
                    rotation_mode = 'anchor'
                    )
    
    # set up the xlimits
    if not keywords['xlimr'] and not keywords['xliml']:
        xl,xr = 2 * [keywords['xlim']]
    else:
        xl,xr = keywords['xliml'],keywords['xlimr']

    # set up the xlimits
    if not keywords['ylimt'] and not keywords['ylimb']:
        yb,yt = 2 * [keywords['ylim']]
    else:
        yb,yt = keywords['ylimb'],keywords['ylimt']

    plt.xlim((min(xvals)-xl,max(xvals)+xr))
    plt.ylim((min(yvals)-yb,max(yvals)+yt))           
    #plt.xlim(min(xvals)-keywords['xlim'],max(xvals)+keywords['xlim'])
    #plt.ylim(min(yvals)-keywords['ylim'],max(yvals)+keywords['ylim'])

    plt.subplots_adjust(
            left = keywords['left'],
            right = keywords['right'],
            top = keywords['top'],
            bottom = keywords['bottom']
            )

    plt.savefig(filename + '.' + fileformat)
    plt.clf()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename+'.'+fileformat))

def plot_concept_evolution(
        scenarios,
        tree,
        fileformat = 'pdf',
        degree = 90,
        verbose = True,
        **keywords
        ):
    """
    Plot the evolution of 
    """
    
    # make defaults
    defaults = dict(
            figsize         = (15,15),
            left            = 0.05,
            top             = 0.95,
            bottom          = 0.05,
            right           = 0.95,
            colormap        = mpl.cm.jet,
            edgewidth       = 5,
            radius          = 2.5,
            outer_radius    = 0.5,
            inner_radius    = 0.25,
            cognates        = '',
            usetex          = False,
            latex_preamble  = False,
            textsize        = 8,
            change          = lambda x:x**1.75,
            xlim            = 0,
            ylim            = 0,
            xlimr           = False,
            xliml           = False,
            ylimt           = False,
            ylimb           = False,
            rootsize        = 10,
            legend          = True,
            legendsize      = 5,
            legendAloc      = 'upper right',
            legendBloc      = 'lower right',
            markeredgewidth = 2.5,
            wedgeedgewidth  = 2,
            gain_linestyle            = 'dotted',
            show_labels     = False,
            loss_linestyle = 'solid',
            ax_linewidth = 0,
            labels = {},
            _prefix = '-   ',
            _suffix = '   -',
            colors = {},
            start = 0,
            filename = 'lingpy-{0}'.format(_timestamp())
            )

    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]
    
    # set filename as variable for convenience
    filename = keywords['filename']

    # XXX customize later XXX
    colormap = keywords['colormap']
   
    # switch backend, depending on whether tex is used or not
    backend = mpl.get_backend()
    if keywords['usetex'] and backend != 'pgf':
        plt.switch_backend('pgf')
    elif not keywords['usetex'] and backend != 'TkAgg':
        plt.switch_backend('TkAgg')

    # check for preamble settings
    if keywords['latex_preamble']:
        mpl.rcParams['pgf.preamble'] = keywords['latex_preamble']

    # make a graph
    graph = nx.Graph()

    # get the tgraph
    tgraph = radial_layout(
            tree,
            degree=degree,
            change=keywords['change'],
            start=keywords['start']
            )

    # get the taxa
    taxa = [n[0] for n in tgraph.nodes(data=True) if n[1]['tip']]

    # set the labels
    labels = {}
    for taxon in taxa:
        if taxon in keywords['labels']:
            labels[taxon] = keywords['labels'][taxon]
        else:
            labels[taxon] = taxon
    
    # get the number of paps in order to get the right colors
    cfunc = np.array(np.linspace(10,256,len(scenarios)),dtype='int')

    if not keywords['colors']:
        colors = dict(
                [
                    (
                        scenarios[i][0],
                        mpl.colors.rgb2hex(
                            colormap(
                                cfunc[i]
                                )
                            )
                        ) for i in range(len(scenarios)
                            )
                        ]
                )
    else:
        colors = keywords['colors']
    
    # get the wedges for the paps
    wedges = {}
    linsp = np.linspace(0,360,len(scenarios)+1)
    for i,scenario in enumerate(scenarios):
        pap = scenario[0]
        theta1,theta2 = linsp[i],linsp[i+1]
        wedges[pap] = (theta1,theta2)
    
    if keywords['legend']:
        
        # set the linestyle for the legend
        if keywords['gain_linestyle'] == 'dotted':
            ls = ':'
        elif keywords['gain_linestyle'] == 'dashed':
            ls = '--'

        legendEntriesA = []
        legendTextA = []
        
        # add stuff for the legend
        for pap,gls in scenarios:
            w = mpl.patches.Wedge(
                    (0,0),
                    1,
                    wedges[pap][0],
                    wedges[pap][1],
                    facecolor = colors[pap],
                    zorder = 1,
                    linewidth=keywords['wedgeedgewidth'],
                    edgecolor='black'
                    )
            legendEntriesA += [w]
            legendTextA += [pap]

        # second legend explains evolution
        legendEntriesB = []
        legendTextB = []
        p = mpl.patches.Wedge(
                (0,0),
                1,
                0,
                360,
                facecolor='0.5',
                linewidth=keywords['wedgeedgewidth'],
                edgecolor='black',
                )
        legendEntriesB += [p]
        legendTextB += ['Loss Event']
        p, = plt.plot(
                0,0,
                ls,
                color='black',
                linewidth=keywords['wedgeedgewidth']
                )
        legendEntriesB += [p]
        legendTextB += ['Gain Event']

        # overwrite stuff
        plt.plot(0,0,'o',markersize=2,zorder=2,color='white')

    # iterate over the paps and append states to the graph
    for pap,gls in scenarios:
        
        # get the graph with the model
        g = gls2gml(
                gls,
                tgraph,
                tree,
                filename = ''
                )

        # iterate over the graph
        for n,d in g.nodes(data=True):
            
            # add the node if necessary
            if n not in graph:
                graph.add_node(n)
            
            # add a pap-dictionary if it's not already there
            if 'pap' not in graph.node[n]:
                graph.node[n]['pap'] = {}

            # add data
            graph.node[n]['pap'][pap] = d['state']
    
    # create the figure
    fig = plt.figure(figsize=keywords['figsize'])
    figsp = fig.add_subplot(111)
    figsp.axes.get_xaxis().set_visible(False)
    figsp.axes.get_yaxis().set_visible(False)

    for s in figsp.spines.values():
        s.set_linewidth(keywords['ax_linewidth'])

    plt.axis('equal')

    xvals = []
    yvals = []

    # iterate over edges first
    for nA,nB in g.edges():
        gA = g.node[nA]['graphics']
        gB = g.node[nB]['graphics']
        xA,yA = gA['x'],gA['y']
        xB,yB = gB['x'],gB['y']

        plt.plot(
                [xA,xB],
                [yA,yB],
                '-',
                color = 'black',
                linewidth=keywords['edgewidth']
                )

    # now iterate over the nodes
    for n,d in graph.nodes(data=True):
        cpaps = d['pap']
        states = list(cpaps.values())
        x,y = g.node[n]['graphics']['x'],g.node[n]['graphics']['y']

        # get z-value which serves as zorder attribute
        try:
            z = 6 * len(self.tree.getConnectingEdges('root',n))
        except:
            z = 0

        xvals += [x]
        yvals += [y]
        
        # plot the default marker
        plt.plot(
                x,
                y,
                'o',
                markersize=keywords['rootsize'],
                color='black',
                zorder=50
                )
        # check for origins in cpaps
        if 'O' in cpaps.values():
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius']+keywords['outer_radius'],
                    0,
                    360,
                    facecolor='white',
                    zorder = 57+z,
                    linewidth= keywords['markeredgewidth'],
                    linestyle = keywords['gain_linestyle'],
                    )
            figsp.add_artist(w)
        elif 'o' in cpaps.values():
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius']+keywords['outer_radius'],
                    0,
                    360,
                    facecolor='white',
                    zorder = 56+z,
                    linewidth=keywords['markeredgewidth'],
                    linestyle='solid',
                    )
            figsp.add_artist(w)
        
        if 'L' in cpaps.values() and 'O' in cpaps.values():
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius']+keywords['outer_radius'],
                    0,
                    360,
                    facecolor='0.5',
                    zorder = 58+z,
                    linewidth = keywords['markeredgewidth'],
                    edgecolor='black',
                    linestyle = keywords['loss_linestyle']
                    )
            figsp.add_artist(w)

        elif "L" in cpaps.values():
            w = mpl.patches.Wedge(
                    (x,y),
                    keywords['radius']+keywords['outer_radius'],
                    0,
                    360,
                    facecolor='0.5',
                    zorder = 59+z,
                    linewidth = keywords['markeredgewidth'],
                    edgecolor='black',
                    )
            figsp.add_artist(w)

        # plot all wedges
        for pap in cpaps:
            
            theta1,theta2 = wedges[pap]
            color = colors[pap]

            # check for characteristics of this pap
            if cpaps[pap] == 'L':

                w = mpl.patches.Wedge(
                        (x,y),
                        keywords['radius'],
                        theta1,
                        theta2,
                        facecolor= color,
                        zorder = 61+z,
                        alpha = 0.25,
                        linewidth = keywords['wedgeedgewidth'],
                        edgecolor='black',
                        linestyle = keywords['loss_linestyle']
                        )
                figsp.add_artist(w)
                
            elif cpaps[pap] == 'o':

                w = mpl.patches.Wedge(
                        (x,y),
                        keywords['radius'],
                        theta1,
                        theta2,
                        facecolor=color,
                        zorder = 61+z,
                        linewidth = keywords['wedgeedgewidth'],
                        edgecolor='black'
                        )
                figsp.add_artist(w)

            elif cpaps[pap] == 'O':

                w = mpl.patches.Wedge(
                        (x,y),
                        keywords['radius'],
                        theta1,
                        theta2,
                        facecolor=color,
                        zorder = 61+z,
                        linewidth = keywords['wedgeedgewidth'],
                        edgecolor='black',
                        linestyle = keywords['gain_linestyle']
                        )
                figsp.add_artist(w)         

        # add the labels if this option is chosen
        if keywords['show_labels']:
            # if node is a tip
            if tgraph.node[n]['tip']:

                # get the values
                gf = tgraph.node[n]['graphics']
                r = gf['angle']
                x,y = gf['x'],gf['y']
                ha = gf['s']
                
                # modify the text
                if ha == 'left':
                    text = keywords['_prefix']+ labels[n]
                else:
                    text = labels[n] + keywords['_suffix']
                
                # plot the text
                plt.text(
                        x,
                        y,
                        text,
                        size = keywords['textsize'],
                        va = 'center',
                        ha = ha,
                        fontweight = 'bold',
                        color = 'black',
                        rotation = r,
                        rotation_mode = 'anchor',
                        zorder = z
                        )
                
    
    # set up the xlimits
    if not keywords['xlimr'] and not keywords['xliml']:
        xl,xr = 2 * [keywords['xlim']]
    else:
        xl,xr = keywords['xliml'],keywords['xlimr']

    # set up the xlimits
    if not keywords['ylimt'] and not keywords['ylimb']:
        yb,yt = 2 * [keywords['ylim']]
    else:
        yb,yt = keywords['ylimb'],keywords['ylimt']

    plt.xlim((min(xvals)-xl,max(xvals)+xr))
    plt.ylim((min(yvals)-yb,max(yvals)+yt))

    prop = mpl.font_manager.FontProperties(size=keywords['legendsize'])
    
    if keywords['legend']:
        legend1 = plt.legend(
                legendEntriesA,
                legendTextA,
                loc=keywords['legendAloc'],
                numpoints=1,
                prop=prop
                )
        plt.legend(
                legendEntriesB,
                legendTextB,
                loc=keywords['legendBloc'],
                prop = prop
                )
        figsp.add_artist(legend1)

    plt.subplots_adjust(
            left= keywords['left'],
            right= keywords['right'],
            top= keywords['top'],
            bottom= keywords['bottom']
            )


    plt.savefig(filename + '.'+fileformat)
    plt.clf()
    if rcParams['verbose']: print(rcParams['M_file_written'].foramt(filename+'.'+fileformat))
