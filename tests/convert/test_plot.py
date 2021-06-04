import pytest

from lingpy.basic.wordlist import Wordlist
from lingpy.convert.plot import *


@pytest.fixture
def Plt(mocker):
    class Plt_(mocker.MagicMock):
        @staticmethod
        def plot(*args, **kw):
            return [mocker.MagicMock()]
    return Plt_()


@pytest.fixture
def Sch(mocker):
    class Sch_(mocker.MagicMock):
        @staticmethod
        def dendrogram(*args, **kw):
            return {'ivl': kw['labels']}
    return Sch_()

@pytest.fixture
def gls():
    return[('Haerbin', 0),
                    ('Beijing', 0),
                    ('Xi’an', 0),
                    ('Lanzhou', 0),
                    ('Kunming', 1),
                    ('Tianjin', 0),
                    ('edge.13', 1),
                    ('edge.11', 0)]


@pytest.fixture
def scenarios(gls):
    return [
        ['p1', gls],
        ['p2', [('Taiyuan', 0),
                    ('Nanning', 0),
                    ('Wenzhou', 0),
                    ('Lanzhou', 0),
                    ('edge.24', 0),
                    ('Pingyao', 0),
                    ('Jinan', 1),
                    ('edge.6', 0),
                    ('Kunming', 0),
                    ('Zhengzhou', 0),
                    ('Xi’an', 0),
                    ('edge.10', 0),
                    ('edge.14', 0),
                    ('Shexian', 0),
                    ('edge.26', 1),
                    ('Suzhou', 0)]]]


@pytest.fixture
def tree():
    return "((((((((Taiyuan,Pingyao,Huhehaote)," \
                    "((((Xi’an,Xining,Zhengzhou)," \
                    "(Lanzhou,Yinchuan,Wulumuqi))," \
                    "(((Tianjin,Jinan),Qingdao),Beijing,Haerbin))," \
                    "(((Guiyang,Kunming),Chengdu,Wuhan)," \
                    "(Nanjing,Hefei)))),(Xiangtan,Changsha)),Nanchang)," \
                    "(Shexian,Tunxi)),((Shanghai,Suzhou,Hangzhou),Wenzhou))," \
                    "(((Xianggang,Guangzhou),Nanning),(Meixian,Taoyuan)))," \
                    "((((Xiamen,Taibei),Shantou,Haikou),Fuzhou),Jian’ou));"


def test_plots(mocker, Plt, Sch, gls, tree, scenarios, tmp_path, test_data):
    mocker.patch('lingpy.convert.plot.mpl', new=mocker.MagicMock())
    mocker.patch('lingpy.convert.plot.plt', new=Plt)
    mocker.patch('lingpy.convert.plot.sch', new=Sch)

    plot_gls(gls, tree, filename=str(tmp_path / 'test'))
    plot_tree(tree, filename=str(tmp_path / 'test'))
    plot_concept_evolution(scenarios, tree, filename=str(tmp_path / 'test'))

    wl = Wordlist(str(test_data /'KSL.qlc'))
    wl.calculate('tree')
    plot_heatmap(wl, filename=str(tmp_path / 'test'), ref="cogid", refB="cogid", steps=1)
