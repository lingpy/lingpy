# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-10-10 16:31
# modified : 2014-02-19 18:19
"""
Basic functions for HTML-plots.
"""

__author__="Johann-Mattis List"
__date__="2014-02-19"


import os
import colorsys
import codecs
import webbrowser
import json

from ..settings import rcParams
from ..read.qlc import read_msa
from ..sequence.sound_classes import pid, token2class

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
        confidence = False,
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
    lingpy.convert.html.msa2html
    lingpy.convert.html.msa2tex

    """
    defaults = dict(
            json = ""
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]
            
    # open the infile
    try:
        data = codecs.open(infile, "r", "utf-8").read()
    except:
        data = codecs.open(infile+'.alm', "r", "utf-8").read()

    # create the outfile
    if not filename:
        filename = rcParams['filename']
    
    # read in the templates
    path = os.path.join(rcParams['_path'],'data','templates')
    
    if main_template:
        html = codecs.open(main_template,'r','utf-8').read()
    else:
        html = codecs.open(
                os.path.join(
                    rcParams['_path'],
                    'data',
                    'templates',
                    'alm2html.html'
                    ),
                'r',
                'utf-8'
                ).read()
    
    if table_template:
        table = codecs.open(table_template,'r','utf-8').read()
    else:
        if not confidence:
            table = codecs.open(
                    os.path.join(
                        rcParams['_path'],
                        'data',
                        'templates',
                        'alm2html.table.html'
                        ),
                    'r',
                    'utf-8'
                    ).read()
        else:
            table = codecs.open(
                    os.path.join(
                        rcParams['_path'],
                        'data',
                        'templates',
                        'alm2html.table.js.html'
                        ),
                    'r',
                    'utf-8'
                    ).read()
    
    # load css
    css = codecs.open(
            os.path.join(
                rcParams['_path'],
                'data',
                'templates',
                'alm.css'
                ),
            'r',
            'utf-8'
            ).read()
    
    # load js
    js = codecs.open(
            os.path.join(
                rcParams['_path'],
                'data',
                'templates',
                'alm.js'
                ),
            'r',
            'utf-8'
            ).read()


    # check for windows-compatibility
    data = data.replace(os.linesep,'\n')[:-1]

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
        confs = []
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
        bas = ' <tr class="{0}{2} taxon" taxon="{3}">\n{1}'

        for tracer,l in enumerate(m):
            # check whether the current line is a borrowing
            if int(l[0]) < 0:
                loan_line = ' loan'
            else:
                loan_line = ''

            # assign the cognate id
            tmp = '  <td>{0}</td>\n'.format(l[0])
            tmp += '  <td>{0}</td>\n'.format(l[1].strip('.'))

            # check alignments for confidence scores
            ipa_string = ''.join([cell.split('/')[0] for cell in
                l[4:]]).replace('-', '')

            tmp += '  <td>{0}</td>\n'.format(ipa_string)
            tmp += '  <td class="{0}">\n'.format(colors[abs(int(l[0]))])
            tmp += '   <table class="{0}">\n'.format(colors[abs(int(l[0]))])
            tmp += '    <tr>\n{0}    </tr>\n   </table>\n  </td>\n </tr>\n'
            
            # check whether another entry follows that is also an alignment,
            # otherwise, there's no need to display a word as an alignment
            cognate_set = False
            if tracer < len(m)-1:
                if abs(int(m[tracer+1][0])) == abs(int(l[0])):
                    cognate_set = True
            if tracer > 0:
                if abs(int(m[tracer-1][0])) == abs(int(l[0])):
                    cognate_set = True
            
            # fill out html for the cognate sets
            if cognate_set:

                alm = ''
                for char in l[4:]:
                    
                    # check for confidence scores
                    if '/' in char:
                        char,conf,num = char.split('/')
                        conf = int(conf) 
                    else:
                        char,conf,rgb = char,(255,255,255),0.0
                    
                    if char == '-':
                        d = 'dolgo_GAP'
                        c = '#bbbbbb'
                    else:
                        d = 'dolgo_'+token2class(char, rcParams['dolgo'])
                        c = token2class(char, rcParams['_color'])

                        # bad check for three classes named differently
                        if d == 'dolgo__':
                            d = 'dolgo_X'
                        elif d == 'dolgo_1':
                            d = 'dolgo_TONE'
                        elif d == 'dolgo_0':
                            d = 'dolgo_ERROR'
                    
                    if confidence:
                        alm += '     '
                        alm += '<td class="char {1}" confidence={0} '.format(
                                conf,
                                d
                                )
                        alm += 'char="{0}" '.format(char)
                        alm += 'onclick="'+"show('{0}')".format(num)+'" '
                        alm += 'num="{0}"'.format(num)
                        alm += '>\n      {0}\n     </td>\n'.format(char)
                    else:
                        alm += '     '
                        alm += '<td class="char {0}">{1}</td>\n'.format(d,char)
            else:
                alm = '      '
                alm += '<td class="{0}">--</td>\n'.format(colors[abs(int(l[0]))])
            
            # format the alignment
            tmp = tmp.format(alm)

            # check for last line, where a new line should be inserted (not the
            # fastest solution, but plotting is not a matter of time, and it
            # suffices it's current purpose
            if tracer < len(m)-1:
                pass
            else:
                if confidence:
                    tmp += ' </table>\n'

                tmp += ' <tr class="empty"><td colspan="4" class="empty">'
                tmp += '<hr class="empty" /></td></tr>\n'
            
            # format the whole string
            tmp_str += bas.format(
                    colors[abs(int(l[0]))],
                    tmp,
                    loan_line,
                    l[1]
                    )

    if not title:
        title = "LexStat - Automatic Cognate Judgments"
    if not shorttitle:
        shorttitle = "LexStat"

    # check for json-attribute
    if keywords['json']:
        keywords['json'] = 'var myjson = '+json.dumps(keywords['json'],
                indent=1)

    html = html.format(
            shorttitle = shorttitle,
            title = title,
            table = tmp_str,
            dataset = dataset,
            javascript = js,
            css = css,
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

def msa2html(
        msa,
        shorttitle = '',
        filename = '',
        template = '',
        **keywords
        ):
    """
    Convert files in ``msa``-format into colored ``html``-format.

    Parameters
    ----------
    msa : dict
        A dictionary object that contains all the information of an MSA object.

    shorttitle : str
        Define the shorttitle of the ``html``-page. If no title is provided,
        the default title ``SCA`` will be used.

    filename : str (default="")
        Define the name of the output file. If no name is defined, the name of
        the input file will be taken as a default.

    template : str (default="")
        The path to the template file. If no name is defined, the basic
        template will be used. The basic template currently used can be found
        under ``lingpy/data/templates/msa2html.html``.

    Examples
    --------
    Load the libary.

    >>> from lingpy import *
    
    Load an ``msq``-file from the test-sets.

    >>> msa = MSA('harry.msq')

    Align the data progressively and carry out a check for swapped sites.

    >>> msa.prog_align()
    >>> msa.swap_check()
    >>> print(msa)
    w    o    l    -    d    e    m    o    r    t
    w    a    l    -    d    e    m    a    r    -
    v    -    l    a    d    i    m    i    r    -

    Save the data to the file ``harry.msa``.

    >>> msa.output('msa',filename='harry')

    Save the ``msa``-object as ``html``.

    >>> msa.output('html',filename='harry')
    
    Notes
    -----
    The coloring of sound segments with respect to the sound class they belong
    to is based on the definitions given in the ``color``
    :py:class:`~lingpy.data.model.Model`. It can easily be changed and adapted.
    

    See also
    --------
    lingpy.convert.html.alm2html
    """
    defaults = dict(
            pid_mode = 1,
            stress = rcParams['stress'],
            css = False,
            css_link = False,
            css_path = 'msa.css'
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]
    
    # while alm-format can be read from the text-file without problems,
    # msa-format should be loaded first (once this is already provided), the
    # loss in speed won't matter much, since output of data is not a daily task
    
    path = os.path.join(rcParams['_path'],'data','templates')

    # load templates
    if not template:
        html = codecs.open(os.path.join(path,'msa2html.html'),'r','utf-8').read()
    else:
        html = codecs.open(template,'r','utf-8').read()

    if not keywords['css']:
        css = codecs.open(os.path.join(path,'msa.css'), 'r', 'utf-8').read()
    else:
        css = codecs.open(keywords['css'], 'r', 'utf-8').read()
    
    # treat the msa-object as a file and try to load the file if this is the
    # case
    if type(msa) == str:
        msa = read_msa(msa, **keywords)
    else:
        raise ValueError('[!] No filename specified.')
    
    # load dataset, etc.
    dataset = msa['dataset']
    
    # calculate pid score, if it is not passed as argument
    if 'pid_score' not in keywords:
        pid_score = 0
        count = 0
        for i,seqA in enumerate(msa['alignment']):
            for j,seqB in enumerate(msa['alignment']):
                if i < j:
                    pid_score += pid(seqA,seqB,mode=keywords['pid_mode'])
                    count += 1
        pid_score = int(100 * pid_score / count+0.5)
    else:
        pid_score = keywords['pid_score']

    infile = msa['infile']
    seq_id = msa['seq_id']

    # define the titles etc.
    if not shorttitle:
        shorttitle = 'SCA'
    
    # determine the length of the longest taxon
    taxl = max([len(t) for t in msa['taxa']])

    # format css file 
    css = css.replace('TAXON_LENGTH', str(taxl * 10))

    out = ''
    tr = '<tr class="msa">{0}</tr>'
    td_taxon = '<td class="taxon">{0}</td>'
    perc = int(80 / len(msa['alignment'][0]) + 0.5)
    td_residue = '<td class="residue {1}">{0}</td>'
    td_swap = '<td class="residue swap {1}">{0}</td>'
    td_unaligned = '<td class="residue noalign {1}">{0}</td>'
    
    # check for swaps in the alignment
    if 'swaps' in msa:
        swaps = []
        for s in msa['swaps']:
            swaps.extend(s)
    else:
        swaps = []

    # check for 
    local = ['*'] * len(msa['alignment'][0])
    if 'local' in msa:
        local = ['.'] * len(msa['alignment'][0])
        for i in msa['local']:
            local[i] = '*'

    # start iteration
    for i,taxon in enumerate(msa['taxa']):
        tmp = ''
        tmp += td_taxon.format(taxon)
        for j,char in enumerate(msa['alignment'][i]):
            if char == '-':
                d = 'dolgo_GAP'
                c = '#bbbbbb'
            else:
                d = 'dolgo_'+token2class(char, rcParams['dolgo'])
                c = token2class(char, rcParams['_color'])

                # bad check for three classes named differently
                if d == 'dolgo__':
                    d = 'dolgo_X'
                elif d == 'dolgo_1':
                    d = 'dolgo_TONE'
                elif d == 'dolgo_0':
                    d = 'dolgo_ERROR'

                        
            if j in swaps:
                tmp += td_swap.format(char, d)
            elif local[j] != '*':
                tmp += td_unaligned.format(char, d)
            else:
                tmp += td_residue.format(char, d)
        out += tr.format(tmp)
    
    if keywords['css_link']:
        css_string = '<link rel="stylesheet" href="{0}" type="text/css">'.format(keywords['css_path'])
    else:
        css_string = '<style>{0}</style>'.format(css)

    html = html.format(
            table = out,
            dataset = dataset,
            pid = pid_score,
            file = infile,
            sequence = seq_id,
            shorttitle = shorttitle,
            width=len(msa['alignment'][0]),
            table_width='{0}'.format(len(msa['alignment'][0])* 50 + 15 * taxl),
            taxa = len(msa['alignment']),
            uniseqs=len(set(msa['seqs'])),
            css = css_string
            )
    
    if not filename:
        filename = rcParams['filename']

    if not filename.endswith('.html'):
        filename = filename+'.html'

    # check, whether the outfile already exists
    outf = codecs.open(filename,'w','utf-8')
    outf.write(html)
    outf.close()
    if rcParams['verbose']: print(rcParams['M_file_written'].format(filename))

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

def msa2tex(
        infile,
        template = '',
        path = '',
        filename = '',
        verbose = True,
        **keywords
        ):
    """
    Convert an MSA to a tabular representation which can easily be used in
    LaTeX documents.
    """
    defaults = dict(
            pid_mode = 1
            )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # while alm-format can be read from the text-file without problems,
    # msa-format should be loaded first (once this is already provided), the
    # loss in speed won't matter much, since output of data is not a daily task
    
    ## get the path to the templates
    path = os.path.join(rcParams['_path'],'data','templates','msa.tex')

    # load msa
    msa = read_msa(infile)

    ## load templates
    if not template:
        tex = codecs.open(path,'r','utf-8').read()
    else:
        tex = codecs.open(template,'r','utf-8').read()

    # calculate pid score, if it is not passed as argument
    if 'pid_score' not in keywords:
        pid_score = 0
        count = 0
        for i,seqA in enumerate(msa['alignment']):
            for j,seqB in enumerate(msa['alignment']):
                if i < j:
                    pid_score += pid(seqA,seqB,mode=keywords['pid_mode'])
                    count += 1
        pid_score = int(100 * pid_score / count+0.5)
    else:
        pid_score = keywords['pid_score']

    dataset = msa['dataset']
    infile = msa['infile']
    seq_id = msa['seq_id']
    
    # determine the length of the longest taxon
    taxl = max([len(t) for t in msa['taxa']])
    
    height = len(msa['alignment'])
    width = len(msa['alignment'][0])
    
    start = r'\tabular{l'+width*'c'+'}\n'
    start += r'\bf\ttfamily Taxon & \multicolumn{'+str(width)+r'}{l}{\bf\ttfamily Alignment}\\'+'\n'
    
    # check for swaps in the alignment
    if 'swaps' in msa:
        swaps = []
        for s in msa['swaps']:
            swaps.extend(s)
    else:
        swaps = []

    body = start
    for i,taxon in enumerate(msa['taxa']):
        body += r'\ttfamily '+taxon.replace('_',r'\_')
        for j,char in enumerate(msa['alignment'][i]):
            if char != '-':
                cls = rcParams['dolgo'][char]
            elif char == '-':
                cls = 'X'
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


def tokens2html(
        string,
        swaps = [],
        tax_len = None,
        path = '',
        template = ''
        ):
    """
    Function converts an (aligned) string into colored html-format.

    """
    # set the tr-line
    tr = '<tr class="msa">\n{0}\n</tr>'

    # get the percentage scaling factor
    perc = int(80 / len(string) + 0.5)

    # get vals for residue and swaps
    td_residue = '<td class="residue" width="50" align="center" bgcolor="{1}">'+\
            '<font color="{2}">{0}</font></td>\n'
    td_swap = '<td class="residue swap" style="border:solid 3px black" width="50"'+\
            'align="center" bgcolor="{1}"><font color="{2}">{0}</font></td>\n'

    # start with filling the taxon
    out = '<table>'
    
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

    return out+'</table>'
