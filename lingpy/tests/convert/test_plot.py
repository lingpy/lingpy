# *-* coding: utf-8 *-*
from __future__ import unicode_literals, print_function, division
from six import text_type
from lingpy.tests.util import WithTempDir 
from lingpy.tests.util import test_data
from lingpy.convert import tree
from lingpy.basic.tree import Tree
from lingpy.convert.plot import *
from lingpy.basic.wordlist import Wordlist

class TestPlot(WithTempDir):

    def setUp(self):

        WithTempDir.setUp(self)
        self.gls =  [('Haerbin', 0),
            ('Beijing', 0),
            ('Xi’an', 0),
            ('Lanzhou', 0),
            ('Kunming', 1),
            ('Tianjin', 0),
            ('edge.13', 1),
            ('edge.11', 0)]

        self.scenarios = [ 
                ['p1', self.gls],
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
        
        self.tree = "((((((((Taiyuan,Pingyao,Huhehaote),((((Xi’an,Xining,Zhengzhou),(Lanzhou,Yinchuan,Wulumuqi)),(((Tianjin,Jinan),Qingdao),Beijing,Haerbin)),(((Guiyang,Kunming),Chengdu,Wuhan),(Nanjing,Hefei)))),(Xiangtan,Changsha)),Nanchang),(Shexian,Tunxi)),((Shanghai,Suzhou,Hangzhou),Wenzhou)),(((Xianggang,Guangzhou),Nanning),(Meixian,Taoyuan))),((((Xiamen,Taibei),Shantou,Haikou),Fuzhou),Jian’ou));"

    
    def test_plot_gls(self):

        plot_gls(self.gls, self.tree, filename=text_type(test_data('test')))


    def test_plot_tree(self):

        plot_tree(self.tree, filename=text_type(test_data('test')))
    
    def test_plot_concept_evolution(self):

        plot_concept_evolution(self.scenarios, self.tree,
                filename=text_type(test_data('test')))

    def test_plot_heatmap(self):
        
        wl = Wordlist(test_data('KSL.qlc'))
        wl.calculate('tree')
        plot_heatmap(wl,
                filename=text_type(test_data('test')), ref="cogid",
                refB="cogid", steps=1)
        plot_heatmap(wl,
                filename=text_type(test_data('test')), ref="cogid",
                normalized='jaccard', steps=1)
        plot_heatmap(wl,
                filename=text_type(test_data('test')), ref="cogid",
                normalized='swadesh', steps=1)


