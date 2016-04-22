"""
Module provides functions for the transformation of text data into visually appealing
format.
"""
from __future__ import unicode_literals, print_function, division

from lingpy.settings import rcParams
from lingpy import log

import numpy as np
import networkx as nx
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except:
    log.missing_module('matplotlib')
    plt = False

try:
    import scipy.cluster.hierarchy as sch
except:
    log.missing_module('scipy')
    sch = False

from lingpy.thirdparty import cogent as cg
from lingpy.convert.tree import nwk2tree_matrix
from lingpy.convert.graph import gls2gml, radial_layout

def plot_gls(
    gls,
    treestring,
    degree=90,
    fileformat='pdf',
    **keywords
):
    """
    Plot a gain-loss scenario for a given reference tree.
    """

    # get kewyords
    defaults = dict(
        figsize=(15, 15),
        left=0.05,
        top=0.95,
        bottom=0.05,
        right=0.95,
        radius=0.5,
        textsize=8,
        edgewidth=5,
        linewidth=2,
        scale_radius=1.2,
        ylim=1,
        xlim=1,
        text=True,
        gain_color='white',
        loss_color='black',
        gain_linestyle='dotted',
        loss_linestyle='solid',
        ax_linewidth=0,
        filename=rcParams['filename']
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

    tgraph = radial_layout(treestring, degree=degree)

    graph = gls2gml(
        gls,
        tgraph,
        tree
    )

    nodes = []

    # assign nodes and edges
    for n, d in graph.nodes(data=True):
        g = d['graphics']
        x = g['x']
        y = g['y']
        s = d['state']

        nodes += [(x, y, s)]

    # now plot the stuff
    fig = plt.figure(figsize=keywords['figsize'])
    figsp = fig.add_subplot(111)
    figsp.axes.get_xaxis().set_visible(False)
    figsp.axes.get_yaxis().set_visible(False)

    # set the axes linewidht
    for s in figsp.spines.values():
        s.set_linewidth(keywords['ax_linewidth'])

    plt.axis('equal')

    for nA, nB in graph.edges():
        xA = graph.node[nA]['graphics']['x']
        xB = graph.node[nB]['graphics']['x']
        yA = graph.node[nA]['graphics']['y']
        yB = graph.node[nB]['graphics']['y']

        plt.plot(
            [xA, xB],
            [yA, yB],
            '-',
            color='black',
            linewidth=keywords['edgewidth'],
            zorder=1
        )

    # now, iterate over nodes
    for x, y, s in nodes:
        if s == 'O':
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'],
                0, 360,
                facecolor=keywords['gain_color'],
                linewidth=keywords['linewidth'],
                linestyle=keywords['gain_linestyle']
            )
        elif s == 'o':
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'] / keywords['scale_radius'],
                0, 360,
                facecolor=keywords['gain_color'],
                linewidth=keywords['linewidth']
            )
        elif s == 'L':
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'],
                0, 360,
                facecolor=keywords['loss_color'],
                linewidth=keywords['linewidth'],
                linestyle=keywords['loss_linestyle']
            )
        else:
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'] / keywords['scale_radius'],
                0, 360,
                facecolor=keywords['loss_color'],
                linewidth=keywords['linewidth']
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
                size=keywords['textsize'],
                color=c,
                va="center",
                ha="center",
                fontweight='bold'
            )

    # set x and y-values
    xvals = [x[0] for x in nodes]
    yvals = [x[1] for x in nodes]

    plt.xlim(min(xvals) - keywords['xlim'], max(xvals) + keywords['xlim'])
    plt.ylim(min(yvals) - keywords['ylim'], max(yvals) + keywords['ylim'])

    plt.subplots_adjust(
        left=keywords['left'],
        right=keywords['right'],
        top=keywords['top'],
        bottom=keywords['bottom']
    )
    plt.savefig(
        filename + '.' + fileformat
    )
    plt.clf()
    log.file_written(filename + '.' + fileformat)


def plot_tree(
    treestring,
    degree=90,
    fileformat='pdf',
    root="root",
    **keywords
):
    """
    Plot a Newick tree to PDF or other graphical formats.

    Parameters
    ----------
    treestring : str
        A string in Newick format.
    degree : int
        Determine the degree of the tree (this determines how "circular" the
        tree will be).
    fileformat : str (default="pdf")
        Select the fileformat to which the tree shall be written.
    filename : str
        Determine the name of the file to which the data shall be written.
        Defaults to a timestamp.
    figsize : tuple (default=(10,10))
        Determine the size of the figure.
    """

    default = dict(
        ax_linewidth=0,
        bg='black',
        bottom=0.05,
        change=lambda x: x ** 1.75,
        edge_list=[],
        figsize=(10, 10),
        filename=rcParams['filename'],
        fontweight='bold',
        frameon=False,
        ha='center',
        labels=[],
        left=0.05,
        linecolor='black',
        linewidth=5,
        no_labels=False,
        node_dict={},
        nodecolor='black',
        nodesize=10,
        right=0.95,
        start=0,
        textcolor='white',
        textsize='10',
        top=0.95,
        usetex=False,
        va='center',
        xlim=5,
        xliml=False,
        xlimr=False,
        ylim=5,
        ylimb=False,
        ylimt=False,
        rotation_mode='anchor',
        latex_preamble=False,
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
        mpl.rcParams['text.latex.unicode'] = True
    elif not keywords['usetex'] and backend != 'TkAgg':
        plt.switch_backend('TkAgg')

    if keywords['latex_preamble']:
        mpl.rcParams['pgf.preamble'] = keywords['latex_preamble']

    # get the tree-graph
    graph = radial_layout(
        treestring,
        degree=degree,
        change=keywords['change'],
        start=keywords['start'],
        root=root
    )

    # create the figure
    fig = plt.figure(figsize=keywords['figsize'])
    figsp = fig.add_subplot(111)
    figsp.axes.get_xaxis().set_visible(False)
    figsp.axes.get_yaxis().set_visible(False)

    for s in figsp.spines.values():
        s.set_linewidth(keywords['ax_linewidth'])

    # plt.axes(frameon=keywords['frameon'])
    plt.axis('equal')
    plt.xticks([])
    plt.yticks([])

    # get xlim and ylim
    xvals, yvals = [], []
    # start iterating over edges
    for nA, nB, d in list(graph.edges(data=True)) + keywords['edge_list']:

        # get the coordinates
        xA = graph.node[nA]['graphics']['x']
        yA = graph.node[nA]['graphics']['y']
        xB = graph.node[nB]['graphics']['x']
        yB = graph.node[nB]['graphics']['y']

        if 'color' in d:
            plt.plot(
                [xA, xB],
                [yA, yB],
                '-',
                **d
            )
        else:
            plt.plot(
                [xA, xB],
                [yA, yB],
                '-',
                color=keywords['linecolor'],
                linewidth=keywords['linewidth'],
            )

    # get the nodes
    for n, d in graph.nodes(data=True):

        g = d['graphics']
        x, y = g['x'], g['y']

        xvals += [x]
        yvals += [y]

        # try to get information from the node-dict
        try:
            settings = {}
            settings.update(keywords['node_dict'][n])
        except:
            settings = {}

        # overwrite the stuff in keywords
        for k in keywords:
            if k not in settings:
                settings[k] = keywords[k]

        if d['label'].startswith('edge') \
                or d['label'].startswith(root) or keywords['no_labels']:
            plt.plot(
                x,
                y,
                'o',
                markersize=settings['nodesize'],
                color=settings['nodecolor'],
                markeredgewidth=settings['linewidth']
            )
        else:
            try:
                label = keywords['labels'][d['label']]
            except:
                label = d['label']
            if 'rotation' in settings:
                r = settings['rotation']
            else:
                r = g['angle']
            plt.text(
                x,
                y,
                label,
                # d['label'],
                color=settings['textcolor'],
                fontweight=settings['fontweight'],
                va=settings['va'],
                ha=g['s'],
                bbox=dict(
                    facecolor=settings['bg'],
                    boxstyle='square,pad=0.2',
                    ec="none",
                ),
                size=settings['textsize'],
                rotation=r,  # g['angle'],
                rotation_mode=settings['rotation_mode']
            )

    # set up the xlimits
    if not keywords['xlimr'] and not keywords['xliml']:
        xl, xr = 2 * [keywords['xlim']]
    else:
        xl, xr = keywords['xliml'], keywords['xlimr']

    # set up the xlimits
    if not keywords['ylimt'] and not keywords['ylimb']:
        yb, yt = 2 * [keywords['ylim']]
    else:
        yb, yt = keywords['ylimb'], keywords['ylimt']

    plt.xlim((min(xvals) - xl, max(xvals) + xr))
    plt.ylim((min(yvals) - yb, max(yvals) + yt))

    plt.subplots_adjust(
        left=keywords['left'],
        right=keywords['right'],
        top=keywords['top'],
        bottom=keywords['bottom']
    )

    plt.savefig(filename + '.' + fileformat)
    plt.clf()
    log.file_written(filename + '.' + fileformat)


def plot_concept_evolution(
    scenarios,
    tree,
    fileformat='pdf',
    degree=90,
    **keywords
):
    """
    Plot the evolution according to the MLN method of all words for a given concept.
    
    Parameters
    ----------
    tree : str
        A tree representation in Newick format.
    fileformat : str (default="pdf")
        A valid fileformat according to Matplotlib.
    degree : int (default=90)
        The degree by which the tree is drawn. 360 yields a circular tree, 180
        yields a tree filling half of the space of a circle.
    """

    # make defaults
    defaults = dict(
        figsize=(15, 15),
        left=0.05,
        top=0.95,
        bottom=0.05,
        right=0.95,
        colormap=mpl.cm.jet,
        edgewidth=5,
        radius=2.5,
        outer_radius=0.5,
        inner_radius=0.25,
        cognates='',
        usetex=False,
        latex_preamble=False,
        textsize=8,
        change=lambda x: x ** 1.75,
        xlim=0,
        ylim=0,
        xlimr=False,
        xliml=False,
        ylimt=False,
        ylimb=False,
        rootsize=10,
        legend=True,
        legendsize=5,
        legendAloc='upper right',
        legendBloc='lower right',
        markeredgewidth=2.5,
        wedgeedgewidth=2,
        gain_linestyle='dotted',
        loss_linestyle='solid',
        ax_linewidth=0,
        labels={},
        _prefix='-   ',
        _suffix='   -',
        colors={},
        start=0,
        filename=rcParams['filename'],
        loss_alpha=0.1,
        loss_background='0.75',
        edges=[],
        hedge_color="black",
        hedge_width=5,
        hedge_linestyle='dashed',
    )
    keywords.update(defaults)

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
    cfunc = np.array(np.linspace(10, 256, len(scenarios)), dtype='int')

    if not keywords['colors']:
        colors = {scenarios[i][0]: mpl.colors.rgb2hex(colormap(cfunc[i]))
                  for i in range(len(scenarios))}
    else:
        colors = keywords['colors']

    # get the wedges for the paps
    wedges = {}
    linsp = np.linspace(0, 360, len(scenarios) + 1)
    for i, scenario in enumerate(scenarios):
        pap = scenario[0]
        theta1, theta2 = linsp[i], linsp[i + 1]
        wedges[pap] = (theta1, theta2)

    if keywords['legend']:

        # set the linestyle for the legend
        if keywords['gain_linestyle'] == 'dotted':
            ls = ':'
        elif keywords['gain_linestyle'] == 'dashed':
            ls = '--'

        legendEntriesA = []
        legendTextA = []

        # add stuff for the legend
        for pap, gls in scenarios:
            w = mpl.patches.Wedge(
                (0, 0),
                1,
                wedges[pap][0],
                wedges[pap][1],
                facecolor=colors[pap],
                zorder=1,
                linewidth=keywords['wedgeedgewidth'],
                edgecolor='black'
            )
            legendEntriesA += [w]
            legendTextA += [pap]

        # second legend explains evolution
        legendEntriesB = []
        legendTextB = []
        p = mpl.patches.Wedge(
            (0, 0),
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
            0, 0,
            ls,
            color='black',
            linewidth=keywords['wedgeedgewidth']
        )
        legendEntriesB += [p]
        legendTextB += ['Gain Event']

        # overwrite stuff
        plt.plot(0, 0, 'o', markersize=2, zorder=2, color='white')

    # iterate over the paps and append states to the graph
    for pap, gls in scenarios:

        # get the graph with the model
        g = gls2gml(
            gls,
            tgraph,
            tree,
            filename=''
        )

        # iterate over the graph
        for n, d in g.nodes(data=True):

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
    for nA, nB in g.edges():
        gA = g.node[nA]['graphics']
        gB = g.node[nB]['graphics']
        xA, yA = gA['x'], gA['y']
        xB, yB = gB['x'], gB['y']

        plt.plot(
            [xA, xB],
            [yA, yB],
            '-',
            color='black',
            linewidth=keywords['edgewidth']
        )

    # add horizontal edges if this option is chosen
    if keywords['edges']:
        # get the coordinates
        for nA, nB in keywords['edges']:
            gA = g.node[nA]['graphics']
            gB = g.node[nB]['graphics']
            xA, yA = gA['x'], gA['y']
            xB, yB = gB['x'], gB['y']

            plt.plot(
                [xA, xB],
                [yA, yB],
                '-',
                color=keywords['hedge_color'],
                linewidth=keywords["hedge_width"],
                linestyle=keywords['hedge_linestyle']
            )

    # now iterate over the nodes
    for n, d in graph.nodes(data=True):
        cpaps = d['pap']
        x, y = g.node[n]['graphics']['x'], g.node[n]['graphics']['y']

        # get z-value which serves as zorder attribute
        try:
            z = 6 * len(tree.getConnectingEdges('root', n))
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
                (x, y),
                keywords['radius'] + keywords['outer_radius'],
                0,
                360,
                facecolor='white',
                zorder=57 + z,
                linewidth=keywords['markeredgewidth'],
                linestyle=keywords['gain_linestyle'],
            )
            figsp.add_artist(w)
        # check for retentions
        elif 'o' in cpaps.values():
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'] + keywords['outer_radius'],
                0,
                360,
                facecolor='white',
                zorder=56 + z,
                linewidth=keywords['markeredgewidth'],
                linestyle='solid',
            )
            figsp.add_artist(w)

        if 'L' in cpaps.values() and 'O' in cpaps.values():
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'] + keywords['outer_radius'],
                0,
                360,
                facecolor=keywords['loss_background'],
                zorder=58 + z,
                linewidth=keywords['markeredgewidth'],
                edgecolor='black',
                linestyle=keywords['loss_linestyle']
            )
            figsp.add_artist(w)

        elif "L" in cpaps.values():
            w = mpl.patches.Wedge(
                (x, y),
                keywords['radius'] + keywords['outer_radius'],
                0,
                360,
                facecolor=keywords['loss_background'],
                zorder=59 + z,
                linewidth=keywords['markeredgewidth'],
                edgecolor='black',
            )
            figsp.add_artist(w)

        # plot all wedges
        for pap in cpaps:

            theta1, theta2 = wedges[pap]
            color = colors[pap]

            # check for characteristics of this pap

            # if it's a loss
            if cpaps[pap] == 'L':

                w = mpl.patches.Wedge(
                    (x, y),
                    keywords['radius'],
                    theta1,
                    theta2,
                    facecolor=color,
                    zorder=61 + z,
                    alpha=keywords['loss_alpha'],  # 0.25,
                    linewidth=keywords['wedgeedgewidth'],
                    edgecolor='black',
                    linestyle=keywords['loss_linestyle']
                )
                figsp.add_artist(w)

            elif cpaps[pap] == 'o':

                w = mpl.patches.Wedge(
                    (x, y),
                    keywords['radius'],
                    theta1,
                    theta2,
                    facecolor=color,
                    zorder=61 + z,
                    linewidth=keywords['wedgeedgewidth'],
                    edgecolor='black'
                )
                figsp.add_artist(w)

            elif cpaps[pap] == 'O':

                w = mpl.patches.Wedge(
                    (x, y),
                    keywords['radius'],
                    theta1,
                    theta2,
                    facecolor=color,
                    zorder=61 + z,
                    linewidth=keywords['wedgeedgewidth'],
                    edgecolor='black',
                    linestyle=keywords['gain_linestyle']
                )
                figsp.add_artist(w)

                # add the labels if this option is chosen
        if keywords['labels']:
            # if node is a tip
            if tgraph.node[n]['tip']:

                # get the values
                gf = tgraph.node[n]['graphics']
                r = gf['angle']
                x, y = gf['x'], gf['y']
                ha = gf['s']

                # modify the text
                if ha == 'left':
                    text = keywords['_prefix'] + labels[n]
                else:
                    text = labels[n] + keywords['_suffix']

                # plot the text
                plt.text(
                    x,
                    y,
                    text,
                    size=keywords['textsize'],
                    va='center',
                    ha=ha,
                    fontweight='bold',
                    color='black',
                    rotation=r,
                    rotation_mode='anchor',
                    zorder=z
                )

    # set up the xlimits
    if not keywords['xlimr'] and not keywords['xliml']:
        xl, xr = 2 * [keywords['xlim']]
    else:
        xl, xr = keywords['xliml'], keywords['xlimr']

    # set up the xlimits
    if not keywords['ylimt'] and not keywords['ylimb']:
        yb, yt = 2 * [keywords['ylim']]
    else:
        yb, yt = keywords['ylimb'], keywords['ylimt']

    plt.xlim((min(xvals) - xl, max(xvals) + xr))
    plt.ylim((min(yvals) - yb, max(yvals) + yt))

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
            prop=prop
        )
        figsp.add_artist(legend1)

    plt.subplots_adjust(
        left=keywords['left'],
        right=keywords['right'],
        top=keywords['top'],
        bottom=keywords['bottom']
    )

    plt.savefig(filename + '.' + fileformat)
    plt.clf()
    log.file_written(filename + '.' + fileformat)


def plot_heatmap(
    wordlist,
    filename="heatmap",
    fileformat="pdf",
    ref='cogid',
    normalized=False,
    refB='',
    **keywords
):
    """
    Create a heatmap-representation of shared cognates for a given wordlist.

    Parameters
    ----------
    wordlist : lingpy.basic.wordlist.Wordlist
        A Wordlist object containing cognate IDs.
    filename : str (default="heatmap")
        Name of the file to which the heatmap will be written.
    fileformat : str (default="pdf")
        A regular matplotlib-fileformat (pdf, png, pgf, svg).
    ref : str (default="cogid')
        The name of the column that contains the cognate identifiers.
    normalized : {bool str} (default=True)
        If set to c{False}, don't normalize the data. Otherwise, select the
        normalization method, choose between:
        
        * "jaccard" for the Jaccard-distance (see :evobib:`Bategelj1995` for
          details), and
        * "swadesh" for traditional lexicostatistical calculation of shared
          cognate percentages.

    cmap : matplotlib.cm (default=matplotlib.cm.jet)
        The color scheme to be used for the heatmap.
    steps : int (default=5)
        The number of steps in which names of taxa will be written to the axes.
    xrotation : int (default=45)
        The rotation of the taxon-names on the x-axis.
    colorbar : bool (default=True)
        Specify, whether a colorbar should be added to the plot.
    figsize : tuple (default=(10,10))
        Specify the size of the figure.
    tree : str (default='')
        A tree passed for the taxa in Newick-format. If no tree is specified,
        the method looks for a tree object in the Wordlist.

    Notes
    -----
    This function plots shared cognate percentages.

    """
    defaults = dict(
        cmap=mpl.cm.jet,
        textsize=5,
        steps=20,
        xrotation=90,
        colorbar=True,
        matrix=False,
        colorbar_label="Shared Cognates",
        figsize=(10, 5),
        colorbar_shrink=0.75,
        colorbar_textsize=10,
        left=0.01,  # rcParams['phybo_xlimr'],
        right=0.95,  # rcParams['phybo_xliml'],
        top=0.95,  # rcParams['phybo_ylimt'],
        bottom=0.01,  # rcParams['phybo_ylimb']
        tree='',
        normalization="jaccard",
        labels={},  # taxon labels passed for the taxa,
        show_tree=True,
        tree_left=0.1,
        tree_bottom=0.1,
        tree_width=0.2,
        height=0.8,
        width=0.8,
        scale=0.075,
        vmax=1.0,
        vmin=0.0
    )
    for k in defaults:
        if k not in keywords:
            keywords[k] = defaults[k]

    # access the reference tree of the wordlist and create a function that
    # orders the taxa accordingly
    if not keywords['tree']:
        try:
            tree = wordlist.tree
        except:
            raise ValueError("[i] No tree could be found")
    else:
        tree = keywords["tree"]

    # check for normalization
    if normalized:
        if normalized not in ["jaccard", "swadesh"]:
            raise ValueError(
                "Keyword 'normalized' must be one of 'jaccard','swadesh',False.")

    # create an empty matrix
    if not normalized:
        matrix = np.zeros((wordlist.width, wordlist.width), dtype=int)
    else:
        matrix = np.zeros((wordlist.width, wordlist.width), dtype=float)

    # create the figure
    fig = plt.figure(figsize=keywords['figsize'])

    # plot the reference tree
    if keywords['show_tree']:
        tree_matrix, taxa = nwk2tree_matrix(tree)
        ax1 = fig.add_axes(
            [
                keywords['left'],
                keywords['bottom'],
                0.25 * keywords['width'],
                keywords['height']
            ]
        )
        # [0.01,0.1,0.2,0.7])
        d = sch.dendrogram(
            np.array(tree_matrix),
            labels=[t for t in taxa],
            orientation='left',

        )
        taxa = d['ivl'][::-1]
        ax1.set_xticks([])
        ax1.set_yticks([])

        left = keywords['left'] + keywords['scale'] * keywords['width']

    else:
        left = keywords['left']
        taxa = tree.taxa

    # start iterating over taxa in order of the reference tree and fill in the
    # matrix with numbers of shared cognates
    if keywords['matrix']:
        matrix = keywords['matrix']
    else:
        for i, taxonA in enumerate(taxa):
            for j, taxonB in enumerate(taxa):
                if i < j:
                    if normalized in [False, "jaccard"]:
                        cogsA = wordlist.get_list(
                            taxa=taxonA,
                            flat=True,
                            entry=ref
                        )
                        cogsB = wordlist.get_list(
                            taxa=taxonB,
                            flat=True,
                            entry=ref
                        )

                        cogsA, cogsB = set(cogsA), set(cogsB)

                        shared = len(cogsA.intersection(cogsB))

                        if normalized:
                            shared = shared / len(cogsA.union(cogsB))
                    else:
                        cogsA = wordlist.get_dict(
                            taxa=taxonA,
                            entry=ref
                        )
                        cogsB = wordlist.get_dict(
                            taxa=taxonB,
                            entry=ref
                        )

                        shared = 0
                        slots = 0

                        # iterate over cognate sets in meaning slots
                        for key in cogsA.keys():
                            # check whether keys are present, we follow the
                            # STARLING procedure in ignoring missing data
                            if key in cogsA and key in cogsB:

                                # check for shared items
                                if [k for k in cogsA[key] if k in cogsB[key]]:
                                    shared += 1
                                slots += 1
                        try:
                            shared = shared / slots
                        except ZeroDivisionError:
                            log.warn(str(
                                [shared, slots, len(cogsA), len(cogsB), taxonA, taxonB]))
                            shared = 0.0

                    matrix[i][j] = shared

                    # if refB is also a possibiltiy
                    if not refB:
                        matrix[j][i] = shared

                elif i > j and refB:
                    if normalized in [False, "jaccard"]:
                        cogsA = wordlist.get_list(
                            taxa=taxonA,
                            flat=True,
                            entry=refB
                        )
                        cogsB = wordlist.get_list(
                            taxa=taxonB,
                            flat=True,
                            entry=refB
                        )

                        cogsA, cogsB = set(cogsA), set(cogsB)

                        shared = len(cogsA.intersection(cogsB))

                        if normalized:
                            shared = shared / len(cogsA.union(cogsB))
                    else:
                        cogsA = wordlist.get_dict(
                            taxa=taxonA,
                            entry=refB
                        )
                        cogsB = wordlist.get_dict(
                            taxa=taxonB,
                            entry=refB
                        )

                        shared = 0
                        slots = 0

                        # iterate over cognate sets in meaning slots
                        for key in cogsA.keys():
                            # check whether keys are present, we follow the
                            # STARLING procedure in ignoring missing data
                            if key in cogsA and key in cogsB:

                                # check for shared items
                                if [k for k in cogsA[key] if k in cogsB[key]]:
                                    shared += 1
                                slots += 1
                        try:
                            shared = shared / slots
                        except ZeroDivisionError:
                            log.warn(str(
                                [shared, slots, len(cogsA), len(cogsB), taxonA, taxonB]))
                            shared = 0.0

                    matrix[i][j] = shared

                elif i == j:
                    cogs = wordlist.get_list(
                        taxa=taxonA,
                        flat=True,
                        entry=ref
                    )
                    if normalized:
                        matrix[i][j] = 1.0
                    else:
                        matrix[i][j] = len(set(cogs))
    ax2 = fig.add_axes(
        [
            left,  # keywords['left']+0.25 * keywords['width']+0.05,
            keywords['bottom'],
            keywords['width'],
            keywords['height']
        ]
    )
    cmap = keywords['cmap'] 

    # [0.15,0.1,0.7,0.7])

    im = ax2.matshow(matrix, aspect='auto', origin='lower',
                     interpolation='nearest', cmap=keywords['cmap'],
                     vmax=keywords['vmax'], vmin=keywords['vmin']
                     )

    # set the xticks
    steps = int(len(taxa) / keywords['steps'] + 0.5)
    start = int(steps / 2 + 0.5)
    idxs = [0] + list(range(start, len(taxa), steps))
    selected_taxa = [taxa[i] for i in idxs]

    # modify taxon names if this is specified
    for i, t in enumerate(selected_taxa):
        if t in keywords['labels']:
            selected_taxa[i] = keywords['labels'][t]

    ax2.set_xticks([])
    ax2.set_yticks([])

    ax1.spines['bottom'].set_color('#ffffff')
    ax1.spines['top'].set_color('#ffffff')
    ax1.spines['left'].set_color('#ffffff')
    ax1.spines['right'].set_color('#ffffff')

    plt.xticks(
        idxs,
        selected_taxa,
        size=keywords['textsize'],
        rotation=keywords['xrotation'],
        rotation_mode="default"
    )
    plt.yticks(
        idxs,
        selected_taxa,
        size=keywords['textsize'],
    )

    if keywords["colorbar"]:
        plt.imshow(matrix, cmap=keywords['cmap'], visible=False, vmax=keywords['vmax'])
        c = plt.colorbar(im, shrink=keywords['colorbar_shrink'])
        c.set_label(keywords["colorbar_label"], size=keywords['colorbar_textsize'])

    plt.subplots_adjust(
        left=keywords['left'],
        right=keywords['right'],
        top=keywords['top'],
        bottom=keywords['bottom']
    )
    plt.savefig(filename + '.' + fileformat)

    f = open(filename + '.matrix', 'w')
    for i, t in enumerate(taxa):
        f.write('{0:20}'.format(t))
        for j, c in enumerate(matrix[i]):
            if not normalized:
                f.write('\t{0:3}'.format(int(c)))
            else:
                f.write('\t{0:.2f}'.format(c))
        f.write('\n')
    f.close()
    log.file_written(filename + '.' + fileformat)
