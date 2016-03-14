# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
from lingpy.cli import parser, pairwise, multiple, settings, lexstat, alignments
from lingpy.tests.util import test_data, WithTempDir
from six import text_type as str

class test_cli(WithTempDir):
    
    def tearDown(self):
        WithTempDir.tearDown(self)

    def test_alignments(self):
        self.setUp()
        args = parser.parse_args(['alignments', '-i',
            test_data('KSL.qlc'), '-c', 'cogid', '--output-file',
            self.tmp_path('alignments').as_posix()])
        alms = alignments(args)
        assert 'German' in alms[1]['taxa']
        assert len(alms[952]['alignment'][0]) == len(alms[952]['alignment'][1])
        
        # test html output and scorer which was computed before
        args = parser.parse_args(['alignments', '-i',
            test_data('KSL3.qlc'), '-c', 'cogid', '--output-file',
            self.tmp_path('alignments').as_posix(), '--format', 'html',
            '--use-logodds'])
        alms = alignments(args)
        assert 'German' in alms[78]['taxa']
        assert len(alms[884]['alignment'][0]) == len(alms[884]['alignment'][1])


        self.tearDown()
    def test_multiple(self):
    
        self.setUp()
        
        # first test, align string, no output, no input
        args = parser.parse_args(['multiple', '-s', 'woldemort', 'waldemar', 'walter'])
        mlt = multiple(args)
        assert len(mlt) == 3
        assert len(mlt[0]) == len('woldemort')
    
        # second test, test output as file, no input, vary method as sca
        args = parser.parse_args(['multiple', '-s', 'woldemort', 'waldemar',
            'walter', '--method', 'sca', '--output-file',
            self.tmp_path('out.msa').as_posix()])
        mlt = multiple(args)
        assert mlt[0] == list('woldemort')
    
        # third test, test output and input
        # second test, test output as file, no input, vary method as sca
        args = parser.parse_args(['multiple', '-i', test_data('harryp.msa'), 
            '--method', 'sca', '--output-file',
            self.tmp_path('out2.msa').as_posix(), '--align-method', 'library'])
        mlt = multiple(args)
        assert len(mlt[0]) == 7
       
        # fourth test, test output and input with method=basic
        args = parser.parse_args(['multiple', '-i', test_data('harryp.msa'), 
            '--method', 'basic', '--output-file',
            self.tmp_path('out2.msa').as_posix()])
        mlt = multiple(args)
        assert len(mlt[0]) == 7
        assert len([x for x in mlt[1][-1] if x != '-']) == 4
        self.tearDown()
    
    def test_pairwise(self):
    
        self.setUp()
    
        # first test, align string, no output, no input
        args = parser.parse_args(['pairwise', '-s', 'woldemort', 'waldemar'])
        pair = pairwise(args)
        assert len(pair) == 3
        assert len(pair[1]) == len('woldemort')
    
        # second test, test output as file, no input, vary method as sca
        args = parser.parse_args(['pairwise', '-s', 'woldemort', 'waldemar', 
            '--method', 'sca', '--output-file',
            self.tmp_path('out.psa').as_posix(), '--distance'])
        pair = pairwise(args)
        assert pair[0] == list('woldemort')
    
        # third test, test output and input
        # second test, test output as file, no input, vary method as sca
        args = parser.parse_args(['pairwise', '-i', test_data('harry_potter.psa'), 
            '--method', 'sca', '--output-file',
            self.tmp_path('out2.psa').as_posix(), '--mode', 'overlap'])
        pair = pairwise(args)
        assert len(pair[0]) == 3
        assert pair[0][0] == list('w-aldemar')
       
        # fourth test, test output and input with method=basic
        args = parser.parse_args(['pairwise', '-i', test_data('harry_potter.psa'), 
            '--method', 'basic', '--output-file',
            self.tmp_path('out2.psa').as_posix()])
        pair = pairwise(args)
        assert len(pair[1][0]) == len('woldemort')+1
        assert [isinstance(x[2], (float,)) for x in pair]
        self.tearDown()
    
    def test_settings(self):
    
        args = parser.parse_args(['settings', '-p', 'lexstat_threshold', 'lexstat_runs',
            'wrgswf'])
        parms = settings(args)
        assert parms[-1][1].startswith('PARAMETER NOT ')
    
    def test_lexstat(self):
        self.setUp()
        args = parser.parse_args(['lexstat', '-i',
            test_data('KSL.qlc'), 
            '--output-file', self.tmp_path('lexstat').as_posix()])
        cogs = lexstat(args)
        assert cogs == 1080
        self.tearDown()
    

