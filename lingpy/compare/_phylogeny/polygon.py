from __future__ import unicode_literals, print_function, division

from lingpy import log

import numpy as np
import networkx as nx
try:
    import matplotlib.patches as mplPatches
except ImportError:
    log.missing_module('matplotlib')

from .convex_hull import convex_hull


def perp(a):
    """
    code for intersection taken from
    http://stackoverflow.com/questions/3252194/numpy-and-line-intersections
    """
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b


# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return
def seg_intersect(nA, nB):
    a1 = np.array(nA[0], dtype='float')
    a2 = np.array(nA[1], dtype='float')
    b1 = np.array(nB[0], dtype='float')
    b2 = np.array(nB[1], dtype='float')

    da = a2 - a1
    db = b2 - b1
    dp = a1 - b1
    dap = perp(da)
    denom = np.dot(dap, db)
    num = np.dot(dap, dp)
    try:
        x, y = (num / denom) * db + b1
    except ZeroDivisionError:
        return False

    # check whether the point is on both lines
    x1 = sorted([a1[0], a2[0]])
    x2 = sorted([b1[0], b2[0]])
    y1 = sorted([a1[1], a2[1]])
    y2 = sorted([b1[1], b2[1]])

    if x1[0] <= x <= x1[1] and x2[0] <= x <= x2[1]:
        if y1[0] <= y <= y1[1] and y2[0] <= y <= y2[1]:
            if [x, y] not in [k.tolist() for k in [a1, a2, b1, b2]]:
                return True
    return False


def getConvexHull(
    nodes,
    color='orange',
    alpha=0.5,
    polygon=True
):
    if len(nodes) < 3:
        return nodes

    points = np.array([[i[0] for i in nodes], [i[1] for i in nodes]])
    ch_points = convex_hull(points, graphic=False)

    if not polygon:
        return ch_points

    return mplPatches.Polygon(ch_points, closed=True, fill=True, color=color, alpha=alpha,
                              lw=0)


def getPolygonFromNodes(
    nodes,
    color='orange',
    alpha=0.5,
):
    # get all lines
    lines = []
    for i, nA in enumerate(nodes):
        for j, nB in enumerate(nodes):
            if i < j:
                lines.append((nA, nB))

    # get all non-intersecting lines
    ni_lines = []
    for i, lineA in enumerate(lines):
        this_line = lineA
        for j, lineB in enumerate(lines):
            if i != j:
                if seg_intersect(lineA, lineB) and lineA[0].tolist() != lineB[0].tolist() and \
                        lineA[1].tolist() != lineB[1].tolist():
                    this_line = False
        if this_line:
            ni_lines.append(this_line)

    # make a graph of the lines
    g = nx.Graph()
    for a, b in ni_lines:
        xA, yA = a
        xB, yB = b
        absX = abs(xA - xB)
        absY = abs(yA - yB)
        if absX == 0:
            d = absY
        elif absY == 0:
            d = absX
        else:
            d = np.sqrt(absX ** 2 + absY ** 2)
        g.add_edge(a, b, weight=d)

    # sort the lines
    start = g.nodes()[0]
    paths = [start]
    while True:
        tmp = g.edge[start].items()
        neighbors = [i[0] for i in tmp if i[0] not in paths]
        if not neighbors:
            break
        distances = [i[1]['weight'] for i in tmp if i[0] not in paths]
        this_neighbor = neighbors[distances.index(max(distances))]
        start = this_neighbor
        paths.append(this_neighbor)

    return mplPatches.Polygon(paths, closed=True, fill=True, color=color, alpha=alpha,
                              lw=0)
