import os

import pytest
import numpy as np

import lingpy
from lingpy.compare._phylogeny.convex_hull import convex_hull
from lingpy.compare._phylogeny.polygon import seg_intersect, getConvexHull, \
    getPolygonFromNodes
from lingpy.compare._phylogeny.utils import get_acs, check_stats, tstats
from lingpy.compare.phylogeny import PhyBo


@pytest.fixture
def Plt(mocker):
    class Plt_(mocker.MagicMock):
        def plot(self, *args, **kw):
            return [mocker.MagicMock()]

        def fill(self, *args, **kw):
            return [mocker.MagicMock()]

        def text(self, *args, **kw):
            return [mocker.MagicMock()]

        def gca(self, *args, **kw):
            return mocker.MagicMock(patches=[])

        def Polygon(self, *args, **kw):
            return [mocker.MagicMock()]
    return Plt_()


@pytest.fixture
def Nx(mocker):
    class Nx_(mocker.MagicMock):
        def Graph(self, *args, **kw):
            return mocker.MagicMock(nodes=lambda **kw: [(mocker.MagicMock(), mocker.MagicMock())])
    return Nx_()


@pytest.fixture
def SPS(mocker):
    class SPS_(mocker.MagicMock):
        mstats = mocker.MagicMock(kruskalwallis=lambda a, b: [1, 1])
    return SPS_()


def test_convex_hull(mocker, Plt):
    try:
        import pylab
    except ImportError:  # pragma: no cover
        return
    mocker.patch('lingpy.compare._phylogeny.convex_hull.p', new=Plt)
    hp = np.array([(2, 1), (3, 1), (2, 10), (5, 6), (10, 1)]).transpose()
    _ = convex_hull(hp, graphic=False)
    _ = convex_hull(hp, graphic=True)


def test_settings():
    assert lingpy.rc('phybo_suffix') == ' -'


def test_seg_intersect():
    line1 = np.array([[0, 0], [1, 0]])
    line2 = np.array([[4, -5], [4, 2]])
    line3 = np.array([[2, 2], [4, 3]])
    line4 = np.array([[6, 0], [3, 4]])
    l12 = seg_intersect(line1, line2)
    l13 = seg_intersect(line3, line4)

    assert not l12
    assert l13


def test_get_convex_hull(mocker, Plt):
    try:
        import matplotlib.patches
    except ImportError:  # pragma: no cover
        return

    mocker.patch('lingpy.compare._phylogeny.polygon.mplPatches', new=Plt)
    hp = np.array([(2, 1), (3, 1), (2, 10), (5, 6), (10, 1)])
    _ = getConvexHull(hp, polygon=False)
    _ = getConvexHull(hp, polygon=True)


def test_get_polygon_from_nodes(mocker, Nx, Plt):
    mocker.patch('lingpy.compare._phylogeny.polygon.nx', new=Nx)

    try:
        import matplotlib.patches
    except ImportError:  # pragma: no cover
        return

    mocker.patch('lingpy.compare._phylogeny.polygon.mplPatches', new=Plt)
    hp = np.array([(-10, -10), (2, 1), (3, 1), (2, 10), (5, 6), (5, 5), (8, 10), (10, 1)])
    getPolygonFromNodes(hp)


@pytest.fixture
def ifile(test_data):
    return str(test_data / 'phybo.qlc')


def test_utils(mocker, SPS, ifile, tmppath):
    try:
        import scipy.stats
    except ImportError:  # pragma: no cover
        return

    mocker.patch('lingpy.compare._phylogeny.utils.sps', new=SPS)
    phy = PhyBo(ifile, output_dir=str(tmppath))
    phy.analyze()
    get_acs(phy, phy.best_model)
    tstats(phy, phy.best_model, return_dists=True)

    check_stats([phy.best_model], phy, filename=str(tmppath / 'test'), pprint=False)
