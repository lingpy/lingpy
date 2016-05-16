import lingpy
from lingpy.compare._phylogeny.convex_hull import area_of_triangle, convex_hull
from lingpy.compare._phylogeny._settings import *
from lingpy.compare._phylogeny.polygon import seg_intersect, getConvexHull, getPolygonFromNodes
from lingpy.compare._phylogeny.utils import get_acs, check_stats, tstats
from lingpy.compare.phylogeny import PhyBo
from lingpy.tests.util import WithTempDir
from lingpy.tests.util import test_data
import os

from mock import MagicMock, patch
import numpy as np

class Plt(MagicMock):
    def plot(self, *args, **kw):
        return [MagicMock()]
    def fill(self, *args, **kw):
        return [MagicMock()]
    def text(self, *args, **kw):
        return [MagicMock()]
    def gca(self, *args, **kw):
        return MagicMock(patches=[])
    def Polygon(self, *args, **kw):
        return [MagicMock()]

class Nx(MagicMock):
    def Graph(self, *args, **kw):
        return MagicMock(nodes=lambda **kw: [(MagicMock(), MagicMock())])

    def generate_gml(self, *args):
        yield ''

class SPS(MagicMock):
    mstats = MagicMock(kruskalwallis=lambda a, b: [1, 1])

class Graph(MagicMock):
    def nodes(self, **kw):
        return [(MagicMock(), dict(label='a', graphics=defaultdict(lambda: 2)))]


@patch('lingpy.compare._phylogeny.convex_hull.p', new=Plt())
def test_convex_hull():
    
    hp = np.array([(2,1), (3,1), (2,10), (5,6), (10,1)]).transpose()
    ch1 = convex_hull(hp, graphic=False)
    ch2 = convex_hull(hp, graphic=True)


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
    
@patch('lingpy.compare._phylogeny.polygon.mplPatches', new=Plt())
def test_get_convex_hull():

    hp = np.array([(2,1), (3,1), (2,10), (5,6), (10,1)])
    gch1 = getConvexHull(hp, polygon=False)
    gch2 = getConvexHull(hp, polygon=True)


@patch('lingpy.compare._phylogeny.polygon.mplPatches', new=Plt())
@patch('lingpy.compare._phylogeny.polygon.nx', new=Nx())
def test_get_polygon_from_nodes():

    hp = np.array([(-10,-10), (2,1), (3,1), (2,10), (5,6), (5,5), (8,10), (10,1)])
    getPolygonFromNodes(hp)


class TestUtils(WithTempDir):
    def setUp(self):
        WithTempDir.setUp(self)
        self.ifile = test_data('phybo.qlc')
    
    @patch('lingpy.compare._phylogeny.utils.sps', new=SPS())
    def test_utils(self):
        phy = PhyBo(self.ifile, output_dir=self.tmp)
        phy.analyze()
        get_acs(phy, phy.best_model)
        tstats(phy, phy.best_model, return_dists=True)

        check_stats([phy.best_model], phy, filename=os.path.join(self.tmp,
            'test'), pprint=False)
